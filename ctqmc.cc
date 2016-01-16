#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>
#include "config.h"
#include "local.h"
#include "det.h"
#include "random.h"
#include "common.h"
#include <omp.h>
//#include "mpi.h"
#include <cstdlib>
#include <cassert>
#include <complex>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

using namespace std;

class Ctqmc{
    //Ctqmc class for Monte Carlo simulation
    Common& common;// common parameters between each classes
    RanGSL& random;// random number generator
    Time_config& t_config;// time configuration
    Time_config t_config_old;// time configuration for storing temporary old configuration
    Local& local;// local trace class
    Det& det;// hybridization determinant class
    //key of map is for perturbation order, the value stores the apears time for each perturbation order
    vector<int> accept;//number of accepted proposal
    vector<int> num_moves;

    public:
    vector<map<int,int > > hist;// histogram of the accumulated perturbation order, index of vector is for flavor,

    // Initialize the Ctqmc class
    Ctqmc(Common& common_, RanGSL& random_, Time_config& t_config_, Local& local_, Det& det_): common(common_), 
    random(random_), t_config(t_config_), local(local_), det(det_), t_config_old(common_.flavor)
    {   clog << "ctqm initialized parameters:" << endl;
        clog << "flavor=" << common.flavor << endl;
        clog << "max_steps=" << common.max_steps << endl;
        clog << "U=" << common.U << endl;
        clog << "ed=" << common.ed << endl;
        clog << "V=" << common.V << endl;
        clog << "beta=" << common.beta << endl;
        clog << "minM=" << common.minM << endl;
        clog << "minD=" << common.minD << endl;
        hist.resize(common.flavor);// resize the histogram to the number of flavor
        accept.resize(common.flavor);
        num_moves.resize(common.flavor);
        for(int i=0; i<common.flavor; i++) {
            accept[i] = 0;// initialize accepted proposal to 0
            num_moves[i] = 0;
        }
    };

    // doing monte carlo sampling
    void sampling();
    // output histogram
    void output_hist();
    // calculate average observables
    void average();

    private:
    // Insert a kink
    void insert_a_kink(int fl_, bool wup_);// first index for selected flavor, second is flag for warmup
    // remove a kink
    void remove_a_kink(int fl_, bool wup_);
};

void Ctqmc::insert_a_kink(int fl_, bool wup_)
{
    int fl = fl_;

    //clog << "copied old config" << endl;
    t_config_old=t_config;

#ifdef DEBUG
    clog << "old config:" << endl;
    t_config_old.print_config(fl);
    clog << endl;
#endif

    // generate random start time
    double ts = common.beta*random();
    assert(ts <= common.beta && ts >= 0.);
    double te;

    pair<bool, int> accept_index;// storing accept condition and insert index for insertion (accept_condition, insert_index)
    pair<double, int> nextts_index;// storing next t start and the index for insertion (next_ts, index_for_next_ts)

    accept_index = t_config.try_insert_start_time(fl, ts);

    while(!accept_index.first) {
        //clog << "ts propose=" << ts << endl;
        ts = common.beta*random();
        assert(ts <= common.beta && ts >= 0.);
        accept_index = t_config.try_insert_start_time(fl, ts);
    }
    //clog << "ts propose=" << setprecision(9) << ts << endl;
    //clog <<"accept?" << accept_index.first << " index=" << accept_index.second << endl;
    t_config.insert_start_time(fl, accept_index.second ,ts);
    nextts_index = t_config.find_next_start_time(fl, accept_index.second);
    //clog <<"next ts=" << nextts_index.first << " next index=" << nextts_index.second << endl;

    // generate random end time in the allowed interval
    double lmax;// = (nextts_index.first -ts);
    if(accept_index.second == -1) {
        te = fmod( ts + common.beta*random(),  common.beta );
        lmax = common.beta;
    }
    else if(nextts_index.second == 0) {
        te = fmod( ts + (nextts_index.first - ts + common.beta)*random(), common.beta );
        lmax = (nextts_index.first - ts + common.beta); 
    }
    else { 
        te = ts + (nextts_index.first -ts)*random();
        lmax = (nextts_index.first -ts); 
    }
#ifdef DEBUG
    clog << "lmax=" << lmax << " te propose=" << te << endl;
    clog << endl;
#endif
    t_config.insert_end_time(fl, accept_index.second, te);

#ifdef DEBUG
    clog << "new config" << endl;
    t_config.print_config(fl);
    clog << endl;
#endif

    // calculate the matrix M and determinant
    double det_ratio = det.calc_insert_det_ratio( fl, common.V ,common.beta , t_config, accept_index.second, accept_index.second );

    // calculate the local trace
    //double local_trace= local.calc_insert_local_trace(fl, t_config, accept_index.second);
    double local_trace= local.calc_local_trace(fl,t_config);

    // Metropolis algorithm
    int k = t_config.get_pertur_order(fl);//perturbation order, already is k+1 in literature 
    double accept_rate = lmax*common.beta/( k )*local_trace*det_ratio;
#ifdef DEBUG
    clog << "det_ratio = "<< det_ratio << "  local_trace = " << local_trace  << "  accept rate = " << accept_rate << endl;
#endif    
    if( abs(accept_rate) >= 1 || random() < abs(accept_rate) ) { // accept new configuration
        //clog << "accept!" << endl;
        det.update_M(fl);
        local.update_trace(fl);
        accept[fl]+=1;//increase accepted proposal
        num_moves[fl] += 1;
        if(hist[fl_].count(k) && !wup_) hist[fl_][k]+=1; //check if the perturb order is in the histogram
        else if(!wup_) hist[fl_][k]=1;// else set histogram to 1
    }
    else { // keep current configuration
        //clog << "reject!" << endl;
        t_config = t_config_old;
        num_moves[fl] += 1;
        if(hist[fl_].count(k-1) && !wup_) hist[fl_][k-1]+=1; //check if the perturb order is in the histogram
        else if(!wup_) hist[fl_][k-1]=1;// else set histogram to 1
    }
    //clog << "Matrix M is" << endl;
    //det.print_M(fl);

};

void Ctqmc::remove_a_kink(int fl_, bool wup_)
{
    int fl = fl_; 
    int pertur_order = t_config.get_pertur_order(fl); 
    int index_to_remove = int( (pertur_order-1)*random());

    //clog << "copied old config" << endl;
    t_config_old=t_config;

#ifdef DEBUG
    clog << "old config:" << endl;
    t_config_old.print_config(fl);
    clog << endl;
#endif

    double lmax;
    double ts_to_remove = t_config.get_t_start_at( fl,index_to_remove );
    double ts_to_remove_next = t_config.get_t_start_at(fl, (index_to_remove+1) % pertur_order );
    if(pertur_order==1) {
        lmax = common.beta;
    }
    else if( ts_to_remove < ts_to_remove_next  ) {
        lmax = ts_to_remove_next - ts_to_remove;
    }
    else {
        lmax = ts_to_remove_next - ts_to_remove + common.beta;
    }
#ifdef DEBUG    
    clog <<"lmax=" <<lmax <<"  remove ts=" << ts_to_remove <<"  and te="<<t_config.get_t_end_at(fl,index_to_remove)<<endl;
#endif
    // remove a kink
    t_config.remove_time_sector(fl, index_to_remove);

#ifdef DEBUG
    clog << "new config:" << endl;
    t_config.print_config(fl);
    clog << endl;
#endif 

    // calculate the matrix M and determinant
    double det_ratio = det.calc_remove_det_ratio( fl, common.V, common.beta , t_config, index_to_remove, index_to_remove );

    // calculate the local trace
    //double local_trace= local.calc_remove_local_trace(fl, t_config_old, index_to_remove);
    double local_trace= local.calc_local_trace(fl,t_config);

    // Metropolis algorithm
    double accept_rate = pertur_order*local_trace*det_ratio/(lmax*common.beta);
#ifdef DEBUG    
    clog << "det_ratio = "<< det_ratio << "  local_trace = " << local_trace  << "  accept rate = " << accept_rate << endl;
#endif
    if( abs(accept_rate) >= 1 || random() < abs(accept_rate) ) { // accept new configuration
        //clog << "accept!" << endl;
        det.update_M(fl);
        local.update_trace(fl);
        accept[fl]+=1;
        num_moves[fl] +=1;
        if(hist[fl_].count(pertur_order-1) && !wup_) hist[fl_][pertur_order-1]+=1; //check if the perturb order is in the histogram
        else if(!wup_) hist[fl_][pertur_order-1]=1;// else set histogram to 1
    }
    else { // keep current configuration
        //clog << "reject!" << endl;
        t_config = t_config_old;
        num_moves[fl] +=1;
        if(hist[fl_].count(pertur_order) && !wup_) hist[fl_][pertur_order]+=1; //check if the perturb order is in the histogram
        else if(!wup_) hist[fl_][pertur_order]=1;// else set histogram to 1
    }

};

void Ctqmc::sampling()
{
    int fl,step;
    double rand;
    bool wup = true;// if warmup or not

    for(int step=0; step<common.max_steps; step++) {
        if (step>common.warmup) wup = false;
#ifndef DEBUG        
        if(step% (common.max_steps/100)==0 && !wup)
#endif            
            clog <<"====================== ctqmc step: " <<step <<"  beta=" << common.beta <<" ========================" << endl;
        rand = random();
        fl = int( random()*common.flavor );//random choose one flavor for insert or delete
         //swithc to the flavor with lower sucess rate       
        if(accept[fl]>accept[(fl+1)%common.flavor]) 
            fl=(fl+1)%common.flavor; 
#ifdef DEBUG        
        clog <<"insert/remove to fl=" << fl << endl;
#endif            

        if(rand > 0.5) {
#ifdef DEBUG            
            clog <<"-------------------------------- insert a kink ------------------------------------" << endl;
#endif       
            insert_a_kink(fl,wup);
//           insert_a_kink(fl,wup);
//           insert_a_kink(fl,wup);
        }
        else {
#ifdef DEBUG       
            clog <<"-------------------------------- remove a kink ------------------------------------" << endl;
#endif           
            if(t_config.get_pertur_order(fl) != 0 )
                remove_a_kink(fl,wup);
        };
    };
    clog <<"========================= MC simulation done! ================================"<< endl;
    for(int i=0; i<common.flavor; i++) {
        clog <<"fl="<<fl <<" accepted move=" << accept[i] << " total MC step=" << num_moves[i] << " accept rate=";
        clog <<double(accept[i])/double(num_moves[i])*100<< "%" <<endl; 
    }
};

void Ctqmc::output_hist(){
    for(int i=0; i<common.flavor; i++) {
        ofstream histout;
        ostringstream convert;
        convert << i;
        string filename = "hist_fl"+convert.str();
        convert.str(" ");
        convert.clear();
        convert << common.beta;
        filename += "_beta"+convert.str()+".dat";
        histout.open(filename.c_str());
        for(map<int,int>::iterator it = hist[i].begin(); it!=hist[i].end(); it++)
           histout << it->first <<"\t" << it->second << endl;
        histout.close();   
    }    
}

void Ctqmc::average(){
    for(int i=0; i<common.flavor; i++) {
        double sum_k=0;
        double sum_steps=0;
        for(map<int,int>::iterator it = hist[i].begin(); it!=hist[i].end(); it++) {
            sum_k += (it->first)*(it->second);
            sum_steps+=it->second;
        }
        clog<<"fl="<<i <<" sum_k="<< sum_k <<" sum_steps="<< sum_steps << " averaged order=" << sum_k/sum_steps<<endl;
    }    
}

//reduce the data from all the threads
void omp_reduce(vector<map<int,int> >& hist_tot_, vector<map<int,int> > hist_, Common& common )
{
    for(int i=0; i<common.flavor; i++) {
        for(map<int,int>::iterator it = hist_[i].begin(); it!=hist_[i].end(); it++) {
           if( !hist_tot_[i].count(it->first) ) hist_tot_[i][it->first]=it->second;
           else hist_tot_[i][it->first]+=it->second;
        }
    }

};

//do statistica averaging
void average(vector<map<int,int> >& hist_, Common& common){
    for(int i=0; i<common.flavor; i++) {
        double sum_k=0;
        double sum_steps=0;
        for(map<int,int>::iterator it = hist_[i].begin(); it!=hist_[i].end(); it++) {
            sum_k += (it->first)*(it->second);
            sum_steps+=it->second;
        }
        clog<<"fl="<<i <<" sum_k="<< sum_k <<" sum_steps="<< sum_steps << " averaged order=" << sum_k/sum_steps<<endl;
    }    
};

void output_hist(vector<map<int,int> >& hist_, Common& common){
    for(int i=0; i<common.flavor; i++) {
        ofstream histout;
        ostringstream convert;
        convert << i;
        string filename = "hist_fl"+convert.str();
        convert.str(" ");
        convert.clear();
        convert << common.beta;
        filename += "_beta"+convert.str()+".dat";
        histout.open(filename.c_str());
        for(map<int,int>::iterator it = hist_[i].begin(); it!=hist_[i].end(); it++)
           histout << it->first <<"\t" << it->second << endl;
        histout.close();   
    }    
}

int main(int argc, char* argv[]) {
    if (argc < 6){
        cerr << "input parameters: max_steps, beta, U, ed, V" << endl;
        return 1;
    }


    int flavor = 2; // flavor(spin)
    int max_steps = atoi(argv[1]);
    double beta = atof(argv[2]); // inverse temperature
    double U = atof(argv[3]); // interaction
    double ed = atof(argv[4]);// impurity level
    double V = atof(argv[5]); // hopping
    double minM = 1e-10;
    double minD = 1e-10;
    int seed = time(0);//1452654131;//time(0);//123456;//time(0);
    int warmup = 50000;
    Common common;
    common.set_params(flavor, max_steps, U, ed, V, beta, minM, minD, warmup);
    vector<map<int,int > > hist_tot(2);

#ifdef _OMP
    #pragma omp parallel private(seed,common)
    {
        clog << "I am thread " << omp_get_thread_num()<<endl;;
        seed = time(0)*5 + 3*omp_get_thread_num();//2321421+ 18*omp_get_thread_num();
#endif        
        // initialize and set the common parameters share between calsses
        // initialize random number generator
        clog << "starting with seed:" << seed << endl;
        RanGSL random(seed); 

        Time_config t_config(common.flavor);// initial time configuration class
        Local local(common.flavor, common.U, common.ed, common.beta);// initialize local trace class
        Det det(common.flavor);// initialize the deteminant class
    
        Ctqmc ctqmc(common, random, t_config, local, det);// initialize ctqmc class
        ctqmc.sampling();// perform Monte Carlo sampling
#ifndef _OMP       
        ctqmc.output_hist();// output histogram
        ctqmc.average();//calculate average observable
#endif        
#ifdef _OMP
        #pragma omp barrier
        #pragma omp critical
        {            
            omp_reduce(hist_tot,ctqmc.hist,common);
        }

    }
    average(hist_tot,common);
    output_hist(hist_tot,common);
#endif

    //MPI_finalize();
    return 0;
}

// declare static parameters in Common class
int Common::flavor;
int Common::max_steps;
double Common::U;    
double Common::ed;
double Common::V;
double Common::beta;
double Common::minM;
double Common::minD;
double Common::warmup;
