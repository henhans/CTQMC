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
    vector<map<int,int > > hist;// histogram of the accumulated perturbation order, index of vector is for flavor,
    //key of map is for perturbation order, the value stores the apears time for each perturbation order
    int accept;//number of accepted proposal

    public:
    // Initialize the Ctqmc class
    Ctqmc(Common& common_, RanGSL& random_, Time_config& t_config_, Local& local_, Det& det_): common(common_), 
    random(random_), t_config(t_config_), local(local_), det(det_), t_config_old(common_.flavor)
    {   /*clog << "ctqm initialized parameters:" << endl;
        clog << "flavor=" << common.flavor << endl;
        clog << "max_steps=" << common.max_steps << endl;
        clog << "U=" << common.U << endl;
        clog << "ed=" << common.ed << endl;
        clog << "t=" << common.t << endl;
        clog << "beta=" << common.beta << endl;
        clog << "minM=" << common.minM << endl;
        clog << "minD=" << common.minD << endl;*/
        hist.resize(common.flavor);// resize the histogram to the number of flavor
        accept = 0;// initialize accepted proposal to 0
    };

    // doing monte carlo sampling
    void sampling();
    // output histogram
    void output_hist();

    private:
    // Insert a kink
    void insert_a_kink(int fl_);
    // remove a kink
    void remove_a_kink(int fl_);
};

void Ctqmc::insert_a_kink(int fl_)
{
    int fl = fl_; //fist consider only one flavor

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
//        clog << "ts propose=" << ts << endl;
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
    double det_ratio = det.calc_insert_det_ratio( fl, common.beta , t_config, accept_index.second, accept_index.second );

    // calculate the local trace
    double local_trace= local.calc_local_trace(t_config);

    // Metropolis algorithm
    int k = t_config.get_pertur_order(fl);//perturbation order, already is k+1 in literature 
    double accept_rate = lmax*common.beta/( k )*local_trace*det_ratio;
#ifdef DEBUG
    clog << "det_ratio = "<< det_ratio << "  local_trace = " << local_trace  << "  accept rate = " << accept_rate << endl;
#endif    
    if( abs(accept_rate) >= 1 || random() < abs(accept_rate) ) { // accept new configuration
        //clog << "accept!" << endl;
        det.update_M(fl);
        local.update_trace();
        accept+=1;//increase accepted proposal
        if(hist[fl_].count(k)) hist[fl_][k]+=1; //check if the perturb order is in the histogram
        else hist[fl_][k]=1;// else set histogram to 1
    }
    else { // keep current configuration
        //clog << "reject!" << endl;
        t_config = t_config_old;
        if(hist[fl_].count(k-1)) hist[fl_][k-1]+=1; //check if the perturb order is in the histogram
        else hist[fl_][k-1]=1;// else set histogram to 1
    }
    //clog << "Matrix M is" << endl;
    //det.print_M(fl);

};

void Ctqmc::remove_a_kink(int fl_)
{
    int fl = fl_; // first consider only one falvor
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
    double det_ratio = det.calc_remove_det_ratio( fl, common.beta , t_config, index_to_remove, index_to_remove );

    // calculate the local trace
    double local_trace= local.calc_local_trace(t_config);

    // Metropolis algorithm
    double accept_rate = pertur_order*local_trace*det_ratio/(lmax*common.beta);
#ifdef DEBUG    
    clog << "det_ratio = "<< det_ratio << "  local_trace = " << local_trace  << "  accept rate = " << accept_rate << endl;
#endif
    if( abs(accept_rate) >= 1 || random() < abs(accept_rate) ) { // accept new configuration
        //clog << "accept!" << endl;
        det.update_M(fl);
        local.update_trace();
        accept+=1;
        if(hist[fl_].count(pertur_order-1)) hist[fl_][pertur_order-1]+=1; //check if the perturb order is in the histogram
        else hist[fl_][pertur_order-1]=1;// else set histogram to 1
    }
    else { // keep current configuration
        //clog << "reject!" << endl;
        t_config = t_config_old;
        if(hist[fl_].count(pertur_order)) hist[fl_][pertur_order]+=1; //check if the perturb order is in the histogram
        else hist[fl_][pertur_order]=1;// else set histogram to 1
    }

};

void Ctqmc::sampling()
{
    int fl = 0;
    for(int step=0; step<common.max_steps; step++) {
#ifndef DEBUG        
        if(step%100==0)
#endif            
            clog <<"====================== ctqmc step: " <<step <<"  beta=" << common.beta <<" ========================" << endl;
        double rand = random();

        if(rand > 0.5) {
            //clog <<"-------------------------------- insert a kink ------------------------------------" << endl;
            insert_a_kink(fl);
//        clog <<"-------------------------------- insert a kink ------------------------------------" << endl;
//           insert_a_kink(fl);
//        clog <<"-------------------------------- insert a kink ------------------------------------" << endl;
//           insert_a_kink(fl);
        }
        else {
            //clog <<"-------------------------------- remove a kink ------------------------------------" << endl;
            if(t_config.get_pertur_order(fl) != 0 )
                remove_a_kink(fl);
        };
    };
    clog <<"========================= MC simulation done! ================================"<< endl;
    clog <<"accepted move=" << accept << " total MC step=" << common.max_steps << " accept rate=";
    clog <<double(accept)/double(common.max_steps)*100<< "%" <<endl; 
};

void Ctqmc::output_hist(){
    for(int i=0; i<common.flavor; i++) {
        ofstream histout;
        ostringstream convert;
        convert << i;
        string filename = "hist"+convert.str()+".dat";
        histout.open(filename.c_str());
        for(map<int,int>::iterator it = hist[i].begin(); it!=hist[i].end(); it++)
           histout << it->first <<"\t" << it->second << endl;
        histout.close();   
    }    
}

int main(int argc, char* argv[]) {
    if (argc < 3){
        cerr << "input parameters: max_steps, beta" << endl;
        return 1;
    }
    int flavor = 1; // flavor(spin)
    int max_steps = atoi(argv[1]);
    double U = 0.0; // interaction
    double ed = -0.0;// impurity level
    double t = 0.5; // hopping
    double beta = atoi(argv[2]); // inverse temperature
    double minM = 1e-10;
    double minD = 1e-10;
    int seed = 1452654131;//time(0);//123456;//time(0);

    // initialize and set the common parameters share between calsses
    Common common;
    common.set_params(flavor, max_steps, U, ed, t, beta, minM, minD);
    // initialize random number generator
    clog << "starting with seed:" << seed << endl;
    RanGSL random(seed); 

    Time_config t_config(common.flavor);// initial time configuration class
    Local local(common.flavor, common.U, common.ed, common.beta);// initialize local trace class
    Det det(common.flavor);// initialize the deteminant class

    Ctqmc ctqmc(common, random, t_config, local, det);// initialize ctqmc class
    ctqmc.sampling();// perform Monte Carlo sampling
    ctqmc.output_hist();// output histogram

    return 0;
}

// declare static parameters in Common class
int Common::flavor;
int Common::max_steps;
double Common::U;    
double Common::ed;
double Common::t;
double Common::beta;
double Common::minM;
double Common::minD;
