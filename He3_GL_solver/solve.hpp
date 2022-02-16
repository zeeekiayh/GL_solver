#ifndef _solve
#define _solve

#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>
#include "ConvergenceAccelerator.hpp"
#include "structures.hpp"
#include "basic_gl_solver.hpp"
// #include "he3bulk.hpp"
#include <chrono> // for timing
#include <complex>
#include <vector>
#include <functional>

using std::cout;
using std::endl;
using std::vector;

using namespace Eigen;
typedef SparseMatrix<dcomplex> SpMat_cd;

// use the relaxation method and Anderson Acceleration to solve
template <class GL_Solver>
VectorXcd Solve_Matrix_Equation(SparseMatrix<dcomplex>&, VectorXcd&, GL_Solver&, in_conditions, vector<int>, string);

#endif