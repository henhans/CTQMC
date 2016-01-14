class Local{
    // Class for local impurity data and operation ( local trace....)
    // local trace is calculated using segment picture
    int flavor;// flavor
    double U; // interaction
    double ed;// impurity level
    double beta; // inverse temperature
    int sign_old; // sign for the last segment -1 if it wrap back to the front, 1 if it doesn't wrap back
    int sign;

    public:
    // initial local class
    Local(int flavor_, double U_, double ed_, double beta_){
        flavor = flavor_;
        U = U_;
        ed = ed_;
        beta = beta_;
        clog << "flavor="<< flavor << " U=" << U <<" ed=" << ed <<" beta=" << beta << endl;
        sign_old = 1;// initialize to 1 cause the empty configuration has sign 1
    };

    double calc_local_trace(Time_config& t_config_);
    void update_trace();
};

double Local::calc_local_trace(Time_config& t_config_){
    int size = t_config_.get_pertur_order(0);
    if(size==0) sign =1;
    else if(t_config_.get_t_start_at(0,size-1) > t_config_.get_t_end_at(0,size-1) ) sign = -1;
    else sign = 1;

    return sign*sign_old;
};

void Local::update_trace() {
    sign_old = sign;
};
