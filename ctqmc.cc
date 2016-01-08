#include <iostream>
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
    Local& local;// local trace class
    Det& det;// hybridization determinant class

    public:
    // Initialize the Ctqmc class
    Ctqmc(Common& common_, RanGSL& random_, Time_config& t_config_, Local& local_, Det& det_): common(common_), 
    random(random_), t_config(t_config_), local(local_), det(det_) 
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
    void insert_a_kink();
    // remove a kink
    void remove_a_kink();
};

void Ctqmc::insert_a_kink()
{
    int fl = 0; //fist consider only one flavor

    // generate random start time
    double ts = common.beta*random();
    assert(ts <= common.beta && ts >= 0.);
    double te;

    pair<bool, int> accept_index;// storing accept condition and insert index for insertion
    pair<double, int> nextts_index;// storing next t start and the index for insertion
    accept_index = t_config.try_insert_start_time(fl, ts);
    while(!accept_index.first) {
        clog << "ts propose=" << ts << endl;
        ts = common.beta*random();
        assert(ts <= common.beta && ts >= 0.);
        accept_index = t_config.try_insert_start_time(fl, ts);
    }
    clog << "ts propose=" << ts << endl;
    clog <<"accept?" << accept_index.first << " index=" << accept_index.second << endl;
    t_config.insert_start_time(fl, accept_index.second ,ts);
    nextts_index = t_config.find_next_start_time(fl, accept_index.second);
    clog <<"next ts=" << nextts_index.first << " next index=" << nextts_index.second << endl;

    // generate random end time in the allowed interval
    if(accept_index.second == -1) te = fmod( ts + common.beta*random(),  common.beta );
    else if(nextts_index.second == 0) te = fmod( ts + (nextts_index.first - ts + common.beta)*random(), common.beta );
    else te = ts + (nextts_index.first -ts)*random();
    clog << "diff=" << (nextts_index.first -ts)<< " te propose=" << te << endl;
    t_config.insert_end_time(fl, accept_index.second, te);

    t_config.print_config(fl);
};

void Ctqmc::remove_a_kink()
{
};

void Ctqmc::sampling()
{
    for(int step=0; step<common.max_steps; step++) {
        clog <<"================================ ctqmc step: " <<step <<" ========================================" << endl;
        insert_a_kink();
    };
};

int main() {
    int flavor = 1; // flavor(spin)
    int max_steps = 3;
    double U = 0.0; // interaction
    double ed = -0.0;// impurity level
    double t = 0.5; // hopping
    double beta = 50; // inverse temperature
    double minM = 1e-10;
    double minD = 1e-10;

    Common common;
    common.set_params(flavor, max_steps, U, ed, t, beta, minM, minD);
    RanGSL random(1234); 

    Time_config t_config(common.flavor);// initial time configuration
    Local local(common.flavor, common.U, common.ed, common.beta);// initial local class
    Det det(common.flavor);

    Ctqmc ctqmc(common, random, t_config, local, det);
    ctqmc.sampling();

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
