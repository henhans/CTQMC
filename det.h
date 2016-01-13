#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
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
    void update_M(int fl_);
};

double Det::calc_insert_det_ratio(int fl_, double beta_, Time_config& t_config_, int row_, int col_) {
    clog << "insert hybridization to row: " << row_ << " col:" << col_ << endl;
    if ( row_ < 0 && col_ <0 ) {
        double te = t_config_.get_t_end_at(fl_,0);
        double ts = t_config_.get_t_start_at(fl_,0);
        Mnew[fl_].resize(1,1);
        if( te > ts) Mnew[fl_](0,0) = hyb_func_flat(beta_, te-ts) ;
        else Mnew[fl_](0,0) = -1.*hyb_func_flat(beta_, te+beta_-ts) ;

#ifdef DEBUG
        clog << "M enlarged is:" << endl;
        clog << Mnew[fl_] << endl;
#endif        
        // return the determinant
        return 1.0/Mnew[fl_](0,0);
    }
    else {
        clog << "old Matrix size = " << M[fl_].rows() <<" config size = "<< t_config_.get_pertur_order(fl_) << endl;
        // ========================================================================================
        // Update the Matrix and calculate the determinant ration using Sherman-Morrison Algorithm
        // ========================================================================================

        // here row->te and col->ts
        MatrixXd D_ri(1, M[fl_].cols()); // Delta[row_,i] row_ is the index for row insertion (fix row)
        MatrixXd D_ic(M[fl_].rows(), 1); // Delta[i,col_] col_ is the index for column insertion(fix col)
        double D_rc; // Delta[row_, col_] the element at the middle 
        MatrixXd L;
        MatrixXd R;
        MatrixXd L2;
        MatrixXd R2;
        double q, p;

        //calculate D_ic, D_ri 
        for( int i=0; i<row_ ; i++) {
            double ts_r = t_config_.get_t_start_at(fl_,i);//ts fix row
            double te_r = t_config_.get_t_end_at(fl_,row_);//te fix row
            double ts_c = t_config_.get_t_start_at(fl_,col_);//ts fix col
            double te_c = t_config_.get_t_end_at(fl_,i);//te fix col
            //clog << "ts_r=" << ts_r << "  te_r=" << te_r << "  ts_c=" << ts_c << "  te_c=" << te_c << endl;
            if ( te_r > ts_r) D_ri(0,i) = hyb_func_flat( beta_, te_r - ts_r );
            else D_ri(0,i) = -1.*hyb_func_flat( beta_, te_r+beta_ - ts_r );
            if ( te_c > ts_c) D_ic(i,0) = hyb_func_flat(beta_, te_c - ts_c);
            else D_ic(i,0) = -1.*hyb_func_flat( beta_, te_c+beta_ - ts_c );
        }
        for( int i=row_+1; i<M[fl_].rows()+1; i++) {
            double ts_r = t_config_.get_t_start_at(fl_,i);//ts fix row
            double te_r = t_config_.get_t_end_at(fl_,row_);//te fix row
            double ts_c = t_config_.get_t_start_at(fl_,col_);//ts fix col
            double te_c = t_config_.get_t_end_at(fl_,i);//te fix col
            //clog << "  ts_r=" << ts_r << "  te_r=" << te_r << "  ts_c=" << ts_c << "  te_c=" << te_c << endl;
            D_ri(0,i-1) = (te_r > ts_r) ? hyb_func_flat( beta_, te_r - ts_r ) : -1.*hyb_func_flat( beta_, te_r+beta_ - ts_r );
            D_ic(i-1,0) = (te_c > ts_c) ? hyb_func_flat( beta_, te_c - ts_c ) : -1.*hyb_func_flat( beta_, te_c+beta_ - ts_c );
            //clog << "D_ri=" << D_ri(0,i-1) << " "<< " D_ic=" << D_ic(i-1,0) << endl;
        }
        double ts = t_config_.get_t_start_at(fl_,col_);
        double te = t_config_.get_t_end_at(fl_,row_);
        D_rc = (te >ts ) ? hyb_func_flat( beta_, te -ts): -1.*hyb_func_flat( beta_, te + beta_ -ts);
#ifdef DEBUG
        clog << "D_ri=" << D_ri << endl;
        clog << "D_ic=" << D_ic << endl;
        clog << "D_rc=" << D_rc << endl;
#else
        clog << "calculating Delta vectors..." << endl;
#endif

        // Constructing require matrix to build new matrix M and calculate the detminant ratio 1/p;
        L = M[fl_]*D_ic;
        R = D_ri*M[fl_];
        q = (D_ri*M[fl_]*D_ic)(0,0);
        p = 1.0 / (D_rc - q);
        L2.resize(L.rows()+1 ,L.cols());
        R2.resize(R.rows() ,R.cols()+1);
        for( int i=0; i<row_; i++) L2(i,0) = L(i,0);
        L2(row_,0) = -1;
        if(row_ < L2.rows()-1) for( int i=row_+1; i<L2.rows(); i++) L2(i,0) = L(i-1,0);
        for( int i=0; i<col_; i++) R2(0,i) = R(0,i);
        R2(0,col_) = -1;
        if(col_ < R2.cols()-1) for( int i=col_+1; i<R2.cols(); i++) R2(0,i) = R(0,i-1);

#ifdef DEBUG
        clog << "L: " << L << endl;
        clog << "R: " << R << endl;
        clog << "p: " << p << endl;      
        clog << "q: " << q << endl;
        clog << "L2: " << L2 << endl;
        clog << "R2: " << R2 << endl;
#else
        clog << "calculating vectors L,R,L2,R2 and numbers p and q..."<<endl;
#endif

        // initialize Mnew[fl_] as zero matrix
        //Mnew[fl_].resize( M[fl_].rows()+1 , M[fl_].cols()+1 );
        Mnew[fl_] = MatrixXd::Zero( M[fl_].rows()+1 , M[fl_].cols()+1 );
        // insert the block using Eigen funcction
        if( row_ == Mnew[fl_].rows()-1 && col_== Mnew[fl_].cols()-1 ) {//inser row and col at the end
            MatrixXd M_u_l = M[fl_].block(0,0,row_,col_);//block matrix at the uper left of M
            //clog <<"M_u_l= " << M_u_l << endl;
            Mnew[fl_].block(0,0,row_,col_) = M_u_l;
        }
        else if( row_ ==0 && col_ ==0 ) {
            MatrixXd M_l_r = M[fl_].block(0,0,row_,col_);//block matrix at the lower right of M
            Mnew[fl_].block(1,1,row_,col_) = M_l_r;
        }
        else {
            MatrixXd M_u_l = M[fl_].block(0,0,row_,col_);//block matrix at the uper left of M
            MatrixXd M_u_r = M[fl_].block(0,col_, row_, M[fl_].cols() - col_ );//block matrix at the uper right of M
            MatrixXd M_l_l = M[fl_].block(row_,0, M[fl_].rows() - row_ ,col_);//block matrix at the lower left of M
            MatrixXd M_l_r = M[fl_].block(row_,col_, M[fl_].rows()-row_, M[fl_].cols()-col_);//block matrix at the lower right of M
            Mnew[fl_].block( 0,0,row_,col_) = M_u_l;
            Mnew[fl_].block( 0,col_+1,row_, M[fl_].cols()-col_ ) = M_u_r;
            Mnew[fl_].block( row_+1,0, M[fl_].rows() - row_ , col_ ) = M_l_l;
            Mnew[fl_].block( row_+1,col_+1, M[fl_].rows() - row_, M[fl_].cols() - col_ ) = M_l_r;
        }

        Mnew[fl_] += p*kroneckerProduct(L2,R2); 
#ifdef DEBUG        
        clog << "M enlarged is: " << endl;
        clog << Mnew[fl_] << endl;
#else
        clog << "constructing Mnew using Sherman-Morrison algorithm" << endl;
#endif

        return 1/p;

    }
};

void Det::update_M(int fl_) {
        M[fl_] = Mnew[fl_];
}
