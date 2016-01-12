#include <iostream>
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
    };

    // doing monte carlo sampling
    void sampling();

    private:
    // Insert a kink
    void insert_a_kink(int fl_);
    // remove a kink
    void remove_a_kink(int fl_);
};

void Ctqmc::insert_a_kink(int fl_)
{
    int fl = fl_; //fist consider only one flavor

    clog << "old config:" << endl;
    t_config_old.print_config(fl);

    // generate random start time
    double ts = common.beta*random();
    assert(ts <= common.beta && ts >= 0.);
    double te;

    pair<bool, int> accept_index;// storing accept condition and insert index for insertion
    pair<double, int> nextts_index;// storing next t start and the index for insertion
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
    clog << "lmax=" << lmax << " te propose=" << te << endl;
    t_config.insert_end_time(fl, accept_index.second, te);

    clog << "new config" << endl;
    t_config.print_config(fl);
    //clog << "copied old config" << endl;
    t_config_old=t_config;
    //t_config_old.print_config(fl);

    // calculate the matrix M and determinant
    double det_ratio;
    det_ratio = det.calc_insert_det_ratio( fl, common.beta , t_config, accept_index.second, accept_index.second );
    clog << "Matrix M is" << endl;
    det.print_M(fl);

    // Metropolis algorithm
    int k = t_config.get_pertur_order(fl);//perturbation order 
    double accept_rate = lmax/( ( k + 1 ) * common.beta )*det_ratio;
    clog << "accept rate = " << accept_rate << endl;
    if( accept_rate >= 1 || random() > accept_rate ) { // accept new configuration
        det.update_M();
    }
    else { // keep current configuration
        t_config = t_config_old;
    }
};

void Ctqmc::remove_a_kink(int fl_)
{
    int fl = fl_; // first consider only one falvor
    int pertur_order = t_config.get_pertur_order(fl);   
    int index_to_remove = int( pertur_order*random());
    
    // remove a kink
    t_config.remove_time_sector(fl, index_to_remove);

    t_config.print_config(fl);
};

void Ctqmc::sampling()
{
    int fl = 0;
    for(int step=0; step<common.max_steps; step++) {
        clog <<"================================ ctqmc step: " <<step <<" ========================================" << endl;
        double rand = random();

//        if(rand > 0.5)
            insert_a_kink(fl);
//        else
//            if(t_config.get_pertur_order(fl) != 0 )
//                remove_a_kink(fl);
    };
};

int main() {
    int flavor = 1; // flavor(spin)
    int max_steps = 2;
    double U = 0.0; // interaction
    double ed = -0.0;// impurity level
    double t = 0.5; // hopping
    double beta = 50; // inverse temperature
    double minM = 1e-10;
    double minD = 1e-10;
    int seed = 1234;//time(0);

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
