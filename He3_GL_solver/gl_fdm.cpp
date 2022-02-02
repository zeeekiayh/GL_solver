#include "solve.hpp"
#include "he3bulk.hpp"
#include "includes.hpp"
#include "complexClasses.hpp"
#include "basic_gl_solver.hpp"
// #include "structures.hpp"
// #include "readwrite.hpp"

// #include <eigen/Eigen/Dense>
// using namespace Eigen;

using namespace std;

int main() {
	const int Nop=5; // number of the order parameter components

	in_conditions cond; // holds all the conditions and parameters
	Bound_Cond eta_BC[Nop]; // boundary conditions for OP components
	complex<double> eta[Nop]; // op vector
	Matrix2d **gradK;          // gradient coefficients in GL functional
	gradK = new Matrix2d *[Nop];
	// initialize sub-matrices in K
	for(int i = 0; i <Nop; i++) gradK[i] = new Matrix2d [Nop];

	// get all the conditions and parameters from the file and confirm
	read_input_data(Nop, cond, eta_BC, gradK,  "conditions.txt");
	confirm_input_data(Nop, cond, eta_BC, gradK);

	VectorXcd guess; // a vector to hold the initial guess
	Three_Component_GL_Solver gls(cond); // the object that builds the solver matrix from K (and BC's)

	// make a sc op object
	FiveCompHe3 fch3(Nop, eta);

	// start timing
	auto start = std::chrono::system_clock::now();

	cout << "Building problem...";
	gls.BuildProblem(gradK, eta_BC); // builds all D matrices and then the solver matrix
	cout << "done" << endl;

	// make a guess
	int size = gls.getSolverMatrixSize(); // the number of columns in the matrix
	guess.resize(size);
	for (int i = 0; i < size; i++) guess(i) = 1.;

	// say how long it has taken so far
	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

	// objects to pass into solver
	vector<int> no_update; // says what elements in rhs (and thus 'guess') should not be changed
	SparseMatrix<dcomplex> SolverMatrix = gls.getSolverMatrix(); //

	// run the solver given the matrix, guess, and sc op object (to get RHS function)
	VectorXcd solution = Solve_Matrix_Equation(SolverMatrix, guess, fch3, cond, no_update, "acceleration");

	// say how long it took to solve
	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

	// save the solution to plot (in python)
	WriteToFile(solution,"data.txt", Nop, cond.SIZEv, cond.SIZEu, cond.STEP);

	/*
	VectorXcd guess_cplx;
	Basic_GL_Solver gls("conditions.txt");

	// start timing
	auto start = std::chrono::system_clock::now();

	cout << "Building problem...";
	gls.BuildProblem();
	cout << "done" << endl;

	// make a guess
	int size = gls.getSolverMatrixSize(); // the number of columns in the matrix
	guess_cplx.resize(size);
	for (int i = 0; i < size; i++) guess_cplx(i) = 1.;

	auto end = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

	gls.Solve(guess_cplx);

	end = std::chrono::system_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

	VectorXcd solution_cplx = gls_cplx.getSolution();
	gls.WriteToFile(solution_cplx,"cplx_data.txt");

	// cout << "Free energy: " << gls.free_energy() << endl;
	*/

	//------  de-allocating gradK array ------------------
	for(int i = 0; i <Nop; i++) delete[] gradK[i]; 
	delete[] gradK; 
	//----------------------------------------------------

	return 0;
}
