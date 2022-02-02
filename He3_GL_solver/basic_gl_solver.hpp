#ifndef _basic_gl_solver
#define _basic_gl_solver
#include "includes.hpp"

// ====================================================
// A class that holds the mesh of OP components, builds
//    the whole problem (matrices, rhs vector), and
//    solves the system using (accelerated) relaxation
class Basic_GL_Solver {
	protected: // for variables and objects that will be used in derived classes
		int mSize; // .................... number of mesh points (size of the D matrices or OP-component vectors)
		int OP_size; // .................. number of OP components
		bool update = 1; // .............. says whether or not to use the no_update vector
		vector<int> no_update; // ........ stores all the indices that will not be modified in the RHS
		string method = "acceleration"; // if Solve will use normal relaxtion or the accelerated
		VectorXcd op_vector; // .......... the OP values along the mesh stored in vector form
		VectorXcd solution; // ........... to store the solution to the GL equ. (in the single vector form)
		SpMat_cd SolverMatrix; // ........ solver matrix
		SpMat_cd Du2, Dw2, Duw; // ....... derivative matrices: ADD D-MATRIX NAMES HERE
		in_conditions cond; // ........... parameters and conditions
		Bound_Cond *eta_BC; // ........... array for OP component boundary condition

	public:
		// CONSTRUCTORS
		Basic_GL_Solver();
		Basic_GL_Solver(in_conditions);

		// has to be virtual for unique definitions in derived classes
		// virtual VectorXcd RHS() = 0;

		// Non-Virtual Functions; to be defined in derived classes
		// dcomplex F_Grad(int, int, int);
		// double free_energy();

		// METHODS
		VectorXcd getSolution() const;
		int getSolverMatrixSize() const;
		SpMat_cd getSolverMatrix() const;
	   
		// ADD CODE IN THIS FUNCTION TO BUILD ADDITIONAL D-MATRICES
		void Build_D_Matrices(); // Build the derivative matrices

		void InsertCoeff_Du2(int, int, int, double, vector<Tr>&); // Insert method for the Du^2 matrix derivatives
		void InsertCoeff_Dw2(int, int, int, double, vector<Tr>&); // insert method for the Dw^2 matrix derivatives
		void InsertCoeff_Duw(int, int, int, double, vector<Tr>&); // insert method for the mixed derivative matrices
		// WRITE ADDITIONAL 'InsertCoeff_D**' METHODS HERE
		
		SparseMatrix<dcomplex> Du2_BD(Bound_Cond, int);
		SparseMatrix<dcomplex> Dw2_BD(Bound_Cond, int); // derivative matrix: 2nd-order of 3rd coordinate (i.e. z)
		SparseMatrix<dcomplex> Duw_BD(Bound_Cond, int); // still needs fixed? ... mixed derivative: of 1st and 3rd coordinates (i.e. x & z)
		// WRITE ADDITIONAL 'D**_BD' METHODS HERE

		SpMat_cd BuildSolverMatrix_Test(Matrix2d**, Bound_Cond*);
		void BuildProblem(Matrix2d**, Bound_Cond*);

		// make an educated guess using the boundary conditions
		// IN PROGRESS....
		VectorXcd makeGuess(VectorXcd&);

		// use the relaxation method and Anderson Acceleration to solve
		// void Solve(VectorXcd& guess);

		void WriteToFile_single_vector(VectorXd&, string);
		void WriteToFile_single_vector(VectorXcd&, string);

}; // GL_solver class
// ====================================================

#endif