#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

//This module is use to calculate hybrdization function for constant bath DOS
struct Par {
    /*parameter passed to integral*/
    double beta;
    double tau;
};

double flat_integrand (double x, void * params) {
     /*integrand for flat band hybridization*/
     Par par = *(Par *) params;
     double beta = par.beta;
     double tau = par.tau;
     double f = exp(x*tau) / ( 1 + exp(x*beta) );
     return f;
}

double hyb_func_flat (double beta_, double tau_) {
      /* hybridization function for flat condunction band*/
      gsl_integration_workspace * w 
      = gsl_integration_workspace_alloc (1000);
            
      double result, error;

      Par par;
      par.beta = beta_;
      par.tau = tau_;

      //printf("%.18f \n",par.beta);
      //printf("%.18f \n",par.tau);

      gsl_function F;
      F.function = &flat_integrand;
      F.params = &par;

      gsl_integration_qags (&F, -1, 1, 0, 1e-7, 1000,
                                  w, &result, &error); 

      //printf ("result          = % .18f\n", result);
      //printf ("estimated error = % .18f\n", error);

      gsl_integration_workspace_free (w);

      return result;
}

/*int main() {
    for (int i=0; i<100; i++)
        printf("%.18f \t %.18f \n",i*50.0/100.,hyb_func_flat(50.0, i*50.0/100.));
}
*/
