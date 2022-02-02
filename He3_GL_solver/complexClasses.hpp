#ifndef _complex_includes
#define _complex_includes

// #include <eigen/Eigen/Dense>
// #include <eigen/Eigen/Sparse>
// #include "structures.hpp"
// #include <complex>
// using std::cout;
// using std::endl;
// using namespace Eigen;
// typedef SparseMatrix<dcomplex> SpMat_cd;
#include "includes.hpp"
#include "basic_gl_solver.hpp"
#include "he3bulk.hpp"

// =======================================================
// Derived class from basic_gl_solver.hpp: Basic_GL_Solver
//  for complex, 3-component GL solver
    class Three_Component_GL_Solver : public Basic_GL_Solver {
        public:
        Three_Component_GL_Solver(in_conditions);
        VectorXcd RHS();
        double free_energy();
    };
// =======================================================

// =======================================================
// Derived class from basic_gl_solver.hpp: Basic_GL_Solver
//  for complex, 5-component GL solver
    class Five_Component_GL_Solver : public Basic_GL_Solver {
        public:
        Five_Component_GL_Solver(in_conditions);
        VectorXcd RHS();
        double free_energy();
    };
// =======================================================

#endif
