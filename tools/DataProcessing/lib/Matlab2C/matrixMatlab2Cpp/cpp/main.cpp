#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>

#include <iostream>
#include <vector>
#include <string>

#include "Matrix.h"

using namespace std;

/** Get test matrix for interpolations. */
template <typename Matrix>
void fillMatrix(Matrix& m) {
    // fill it (source wiki bicubic interp.)
    m(0,0) = 1; m(1,0) = 2; m(2,0) = 4; m(3,0) = 1; m(4,0) = 1;
    m(0,1) = 6; m(1,1) = 3; m(2,1) = 5; m(3,1) = 2; m(4,1) = 4;
    m(0,2) = 4; m(1,2) = 2; m(2,2) = 1; m(3,2) = 5; m(4,2) = 4;
    m(0,3) = 5; m(1,3) = 4; m(2,3) = 2; m(3,3) = 3; m(4,3) = 4;
    m(0,4) = 3; m(1,4) = 3; m(2,4) = 3; m(3,4) = 2; m(4,4) = 3;
}

int main(int argc, char * const argv[]) {
    
	typedef Matrix::D2D<double> MatrixD2D;
    
    // generate it
    MatrixD2D m(5,5);
    fillMatrix(m);
    
    // dump it to file & read
    m.dump("matrix.dat");
    MatrixD2D m2("matrix.dat");
    m.load("matrix.dat");
    
    // fancy dump (float and int)
    m.dump("matrixF.dat", Matrix::TID_FLOAT);
	m.dump("matrixI.dat", Matrix::TID_INT);
    
    // show matrix
    for (int j = 0; j < m.getSizeY(); ++j) {
        for (int i = 0; i < m.getSizeX(); ++i) {
            cout << m(i,j);
            if (i < m.getSizeX()-1) cout << ", ";
        }
        cout << endl;
    }
    
    
    return 0;
}
