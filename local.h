class Local{
    // Class for local impurity data and operation ( local trace....)
    // local trace is calculated using segment picture. Currently only the Anderson model is implemented
    int flavor;// flavor
    double U; // interaction
    double ed;// impurity level
    double beta; // inverse temperature
    vector<int> sign_old; // sign for the last segment -1 if it wrap back to the front, 1 if it doesn't wrap back
    vector<int> sign; // each flavor has it's own sign
    //vector<double> l_seg;// the accumulated length for the sectors in each flavor
    //double l_ov;// the loverlap between two flavors

    public:
    // initial local class
    Local(int flavor_, double U_, double ed_, double beta_){
        flavor = flavor_;
        U = U_;
        ed = ed_;
        beta = beta_;
        clog << "flavor="<< flavor << " U=" << U <<" ed=" << ed <<" beta="<< beta << endl;
        sign_old.resize(flavor);
        sign.resize(flavor);
        sign.resize(flavor);
        //l_seg.resize(flavor);
        for(int i=0; i<flavor; i++) sign_old[i] = 1;// initialize to 1 cause the empty configuration has sign 1
    };

    double calc_insert_local_trace(int fl_, Time_config& t_config_, int index_);
    double calc_remove_local_trace(int fl_, Time_config& t_config_, int index_);
    double calc_local_trace(int fl_, Time_config& t_config_);

    void update_trace(int fl_);
};

double Local::calc_insert_local_trace(int fl_, Time_config& t_config_, int index_){
    int fl_other = (fl_+1)%flavor;
    int size = t_config_.get_pertur_order(fl_);
    int size_other = t_config_.get_pertur_order(fl_other);
    double exp_ed, exp_U, lseg, lov;
    double ts = t_config_.get_t_start_at(fl_,index_);
    double te = t_config_.get_t_end_at(fl_,index_);

    if(size==0) { 
        lseg=0;
        sign[fl_] =1;
    }
    else if(t_config_.get_t_start_at(fl_,size-1) > t_config_.get_t_end_at(fl_,size-1) ){
        lseg = (te>ts) ? te-ts : te - ts + beta;
        sign[fl_] = -1;
    }
    else {
        lseg = (te>ts) ? te-ts : te - ts + beta;
        sign[fl_] = 1;
    }

    if(size_other==0) lov=0;
    else{
        int idx=-1;
        bool found=false;
        while( !found && idx < size_other) {
            idx++;
            if( t_config_.get_t_end_at(fl_other,idx) > ts ) found = true;
            //clog << "idx=" << idx<<endl;
        }
        //clog << "found=" <<found << " idx=" << idx<< endl;
        if(t_config_.get_t_end_at(fl_other,idx) < te && t_config_.get_t_start_at(fl_other,idx)< ts)
            lov = t_config_.get_t_end_at(fl_other,idx) - ts;
        else if(t_config_.get_t_end_at(fl_other,idx) < te && t_config_.get_t_start_at(fl_other,idx)> ts)
            lov = t_config_.get_t_end_at(fl_other,idx) - t_config_.get_t_start_at(fl_other,idx);
        else if(t_config_.get_t_end_at(fl_other,idx) > te && t_config_.get_t_start_at(fl_other,idx)< ts)
            lov = te - ts; 
        else if(t_config_.get_t_end_at(fl_other,idx) > te && t_config_.get_t_start_at(fl_other,idx)> ts)
            lov = te - t_config_.get_t_start_at(fl_other,idx);

        if(lov<0) lov=0;
    }


#ifdef DEBUG
    clog << "insert at index="<< index_<<" ts="<< ts<< " te="<< te << endl;
    t_config_.print_config(fl_);
    t_config_.print_config(fl_other);
    clog << "lseg=" << lseg << endl;
    clog << "lov=" << lov << endl;
#endif

    exp_ed = exp(-ed*lseg);
    exp_U = exp(-U*lov);

    return sign[fl_]*sign_old[fl_]*exp_ed*exp_U;
};

double Local::calc_remove_local_trace(int fl_, Time_config& t_config_, int index_){
    int fl_other = (fl_+1)%flavor;
    int size = t_config_.get_pertur_order(fl_);
    int size_other = t_config_.get_pertur_order(fl_other);
    double exp_ed, exp_U, lseg, lov;
    double ts = t_config_.get_t_start_at(fl_,index_);
    double te = t_config_.get_t_end_at(fl_,index_);

    if(size==1) {
        lseg=0;
        sign[fl_] =1;
    }
    else if(t_config_.get_t_start_at(fl_,size-1) > t_config_.get_t_end_at(fl_,size-1) ){
        lseg = (te>ts) ? te-ts : te - ts + beta;
        sign[fl_] = -1;
    }
    else {
        lseg = (te>ts) ? te-ts : te - ts + beta;
        sign[fl_] = 1;
    }

    if(size_other==0) lov=0;
    else{
        int idx=-1;
        bool found=false;
        while( !found && idx < size_other) {
            idx++;
            if( t_config_.get_t_end_at(fl_other,idx) > ts ) found = true;
            //clog << "idx=" << idx<<endl;
        }
        //clog << "found=" <<found << " idx=" << idx<< endl;
        if(t_config_.get_t_end_at(fl_other,idx) < te && t_config_.get_t_start_at(fl_other,idx)< ts)
            lov = t_config_.get_t_end_at(fl_other,idx) - ts;
        else if(t_config_.get_t_end_at(fl_other,idx) < te && t_config_.get_t_start_at(fl_other,idx)> ts)
            lov = t_config_.get_t_end_at(fl_other,idx) - t_config_.get_t_start_at(fl_other,idx);
        else if(t_config_.get_t_end_at(fl_other,idx) > te && t_config_.get_t_start_at(fl_other,idx)< ts)
            lov = te - ts; 
        else if(t_config_.get_t_end_at(fl_other,idx) > te && t_config_.get_t_start_at(fl_other,idx)> ts)
            lov = te - t_config_.get_t_start_at(fl_other,idx);

        if(lov<0) lov=0;
    }


#ifdef DEBUG
    clog << "lseg=" << lseg << endl;
    clog << "lov=" << lov << endl;
#endif

    exp_ed = exp(ed*lseg);
    exp_U = exp(U*lov);

    return sign[fl_]*sign_old[fl_]*exp_ed*exp_U;
};


double Local::calc_local_trace(int fl_, Time_config& t_config_){
    int size = t_config_.get_pertur_order(fl_);

    if(size==0) { 
        sign[fl_] =1;
    }
    else if(t_config_.get_t_start_at(fl_,size-1) > t_config_.get_t_end_at(fl_,size-1) ){
        sign[fl_] = -1;
    }
    else {
        sign[fl_] = 1;
    }

    return sign[fl_]*sign_old[fl_];
};

void Local::update_trace(int fl_) {
    sign_old[fl_] = sign[fl_];
};
