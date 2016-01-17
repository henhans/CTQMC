# CTQMC for practice
Continuous time monte carlo(CTQMC) hybridization algorithm for Anderson model segment update. The only implemented observable is perturbation order. The output is the perturbation order histogram for each spin. At low temperature beta>500 the perturbation order historgram is distorted...probably need to add more update methods(global, anti-segment)...debugging. 

Data structure:

    config.h: Class for storing and updating time configurations. Currently only add one kink and remove one kink is implemented.
    det.h: Class for storing and updating local matrices using Sherman-Morrison Algorithm.
    local.h: Class for calculating local trace by segment picture.

Libraries:

    1.Eigen
    2.gsl
    3.openmp

