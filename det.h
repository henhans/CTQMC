#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Det {
    // class for hybridization matrix delta and its inverse M, and the 
    // routine that calculates the determinant.
    int flavor;
    vector<Matrix2d> M;//Matrix M for inverse hybnridization matrix

    public:
    // initialize the determinant
    Det(int flavor_) {
        flavor = flavor_;
        for( int i=0; i<flavor; i++){
            clog << "setup matrix M for flavor:" << i << endl;
            Matrix2d tmp;
            M.push_back(tmp);
        };
    };

};
