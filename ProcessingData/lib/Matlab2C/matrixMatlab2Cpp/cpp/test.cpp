#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>

#include <iostream>
#include <vector>
#include <string>

#include "Matrix.h"

using namespace std;

/* ==================
 AIM: 
1) read matrix in binary format stored in matlab
2) read and store it in c++
3) print out for check
======================= */

int main(int argc, char * const argv[]) {

	typedef Matrix::D3D<double> MatrixD3D;
    
    /* 1) load matrix from file */
    MatrixD3D m("../testdata/low_anatomy.dat");
	//m.load("../tesdata/test_matrix.dat");
	
	cout<<"size"<<m.getSizeX()<<"x"<<m.getSizeY()<<"x"<<m.getSizeZ()<<endl;
	double a = 0.0;
	double x = 1./a;
	double y = 1./x;
	
	cout << "x="<<x<<" y= "<<y << endl;
    // show matrix
//    for (int j = 0; j < m.getSizeY(); ++j) {
//        for (int i = 0; i < m.getSizeX(); ++i) {
//            cout << m(i,j);
//            if (i < m.getSizeX()-1) cout << ", ";
//        }
//        cout << endl;
//    }
	
	return 0;
}