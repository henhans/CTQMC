#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>
#include "config.h"
#include <omp.h>
#include <cstdlib>
#include <cassert>
#include <complex>

int flavor = 2; // flavor(spin)
double U = 0.0; // interaction
double ed = -0.5;// impurity level
double t = 0.5; // hopping
double beta = 50; // inverse temperature

/*class Ctqmc{
  ctqmc(){
    
  };
};*/

int main() {

  Time_config t_config(flavor);// seup time configuration
  Ctqmc ctqmc(flavor, U, t, beta)

  return 0;
}
