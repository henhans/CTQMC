#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>
#include "config.h"
#include "local.h"
#include <omp.h>
#include <cstdlib>
#include <cassert>
#include <complex>
#include <queue>

using namespace std;

/*class Ctqmc{
  ctqmc(){
    
  };
};*/

int main() {
    int flavor = 2; // flavor(spin)
    double U = 0.0; // interaction
    double ed = -0.5;// impurity level
    double t = 0.5; // hopping
    double beta = 50; // inverse temperature

    Time_config t_config(flavor);// initial time configuration
    Local local(flavor, U, ed, beta);// initial local class
    //Ctqmc ctqmc(flavor, U, t, beta)

    return 0;
}
