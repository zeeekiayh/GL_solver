#ifndef _includes
#define _includes

#include <iostream>
#include <fstream> // for file in/out
#include <string>
#include <math.h>
#include <complex>
#include <vector>
#include <chrono> // for timing
#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>
#include <eigen/unsupported/Eigen/KroneckerProduct> // for the Kronecker product
#include "structures.hpp"
#include "readwrite.hpp"
// #include "he3bulk.hpp"

using std::cout;
using std::endl;
using std::cin;
using std::complex;
using std::ostream;
using std::vector;
using std::string;
using namespace Eigen;

typedef Triplet<double> Tr;
typedef SparseMatrix<dcomplex> SpMat_cd;

// returns the unique id corresponding to each op-component in the mesh.
int ID(int, int, int, int, int);

// Insert the matrix "spMat" into the "location" (i,j) in
//   a sparse matrix of size "size" and return the matrix
SpMat_cd Place_subMatrix(int, int, int, SpMat_cd);

// A special function to access OP elements from a vector that
//   has been flattened from a matrix of size 'sizeU' X 'sizeV'
VectorXcd matrix_operator(VectorXcd& vec, int v, int u, int sizeV, int sizeU, int OPsize, bool debug = false);

// Does these need to be here??
// void Solve_Matrix_Equation(SpMat_cd&, VectorXcd&, VectorXcd);
// #include "complexClasses.hpp" // Now include the real and complexs class header files, since the parent class are defined

#endif
