 class Common {
     //class for common value to pass to all the other class
     public:
     static int flavor;
     static int max_steps;
     static double U;    
     static double ed;
     static double V;
     static double beta;
     static double minM;// The smallest Ddeterminant to accept in order to prevent singularity.
     static double minD;// The smallest Ddeterminant to accept in order to prevent singularity.
     static double warmup;

     static void set_params(int flavor_, int max_steps_ ,double U_, double ed_, double V_, double beta_, double minM_,
     double minD_, int warmup_) {
         flavor = flavor_;
         max_steps = max_steps_;
         U = U_;
         ed = ed_;
         V = V_;
         beta = beta_;
         minM = minM_;
         minD = minD_;
         warmup = warmup_;
     };
 };

