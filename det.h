#include <eigen3/Eigen/Dense>
#include "hybfunc.cc"

using namespace Eigen;

class Det {
    // class for hybridization matrix delta and its inverse M, and the 
    // routine that calculates the determinant, update the matrix...
    int flavor;
    vector<double> det;//
    vector<MatrixXd> M;// Matrix M for inverse hybnridization matrix for each flovrs
    vector<MatrixXd> Mnew;// storing matrix M in previous mc step for each flavors

    public:
    // initialize the determinant
    Det(int flavor_) {
        flavor = flavor_;
        for( int i=0; i<flavor; i++){
            clog << "setup matrix M for flavor:" << i << endl;
            MatrixXd tmp;
            M.push_back(tmp);
            Mnew.push_back(tmp);
        };
    };

    void print_M(int fl_) { clog << M[fl_] << endl;};

    // update determinant by insert a kink(for now row=col)
    double calc_insert_det_ratio(int fl_, double beta_, Time_config& t_config, int row_, int col_);

    // update determinat by remove a kink(for now row=col)
    double calc_remove_det_ratio(int fl_, double beta_, Time_config& t_config, int col_);

    // update matrix M when the move is accept
    void update_M();
};

double Det::calc_insert_det_ratio(int fl_, double beta_, Time_config& t_config_, int row_, int col_) {
    clog << "insert hybridization to row: " << row_ << " col:" << col_ << endl;
    if ( row_ < 0 && col_ <0 ) {
        double te = t_config_.get_t_end_at(fl_,0);
        double ts = t_config_.get_t_start_at(fl_,0);
        M[fl_].resize(1,1);
        if( te > ts) M[fl_](0,0) = hyb_func_flat(beta_, te-ts) ;
        else M[fl_](0,0) = hyb_func_flat(beta_, te+beta_-ts) ;

        // return the determinant
        return 1.0/M[fl_](0,0);
    }
    else {
        clog << "old Matrix size = " << M[fl_].rows() <<" config size = "<< t_config_.get_pertur_order(fl_) << endl;
        MatrixXd D_ri(M[fl_].rows(),1); // Delta[row_,i] row_ is the index for row insertion
        MatrixXd D_ic(1,M[fl_].cols()); // Delta[i,col_] col_ is the index for column insertion
        double D_rc; // Delta[row_, col_] the element at the middle 
        MatrixXd L;
        MatrixXd R;
        MatrixXd L2;
        MatrixXd R2;
        double q, p;

        //calculate D_ic, D_ri 
        for( int i=0; i<row_ ; i++) {
            double ts = ;
            double te = ;
            D_ri(i,1) = hyb_func_flat(beta_,);
            D_ic(1,i) = hyb_func_flat(beta_,);
        }
        for( int i=row_+1; i<M[fl_].rows(); i++) {

        }

        return 1/p;
    }
};

void Det::update_M() {

}
