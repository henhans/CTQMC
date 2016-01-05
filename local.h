class Local{
    // Class for local impurity data and operation ( local trace....)
    // local trace is calculated using segment picture
    int flavor;// flavor
    double U; // interaction
    double ed;// impurity level
    double beta; // inverse temperature
     
    public:
    // initial local class
    Local(int flavor_, double U_, double ed_, double beta_){
        flavor = flavor_;
        U = U_;
        ed = ed_;
        beta = beta_;
        clog << "flavor="<< flavor << " U=" << U <<" ed=" << ed <<" beta=" << beta << endl;
    };
};
