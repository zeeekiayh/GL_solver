// copied Solve() function from Izek's  includes.hpp
// and deleted bunch of lines 
// ------------ maybe not working anymore -----------------

#include <iostream>
#include <fstream> // for file in/out
#include <string>
#include <math.h>
#include <complex>
#include <vector>
#include <chrono> // for timing
#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>
#include "ConvergenceAccelerator.hpp"
#include "structures.hpp"

#include "SC_classes.hpp"

using namespace Eigen;
using namespace std;

// solves self-consistency for linear system of equations
//     M*f = RHS(f)
// RHS function is given through *SC class
// on input f is initial guess, 
// on output it is the required solution 
//---------------------------------------------------------------------------------
void Solver(VectorXcd & f, SpMat_cd M, VectorXcd rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC)
{
	int grid_size=cond.SIZEu * cond.SIZEv; 
	int vect_size=cond.Nop * grid_size; 

	if (no_update.size()) cout << "using no_update" << endl;

	// use LU decomposition to solve the system
	SparseLU< SpMat_cd, COLAMDOrdering<int> > solver;
   	solver.analyzePattern( M ); 
	solver.factorize( M );

	// check to see if the solver failed
	if (solver.info() == Eigen::Success) cout << "\tSolver: successfully built" << endl;
	else if (solver.info() == Eigen::NumericalIssue) { // for debugging non-invertable matrices
		cout << "Solver: numerical issues" << endl;
		// for (int k=0; k < SolverMatrix.outerSize(); ++k) for (SparseMatrix<Scalar_type>::InnerIterator it(SolverMatrix,k); it; ++it) cout << "(" << it.row() << "," << it.col() << ")\t";
		return;
	}

	int cts = 0; // count loops
	// cout << "here ... 1" << endl;
	double err;  // to store current error
	// cout << "here ... 2" << endl;
	VectorXd dummy(vect_size); // dummy free energy variable for RHS function
	// cout << "here ... 3" << endl;
	VectorXcd df(vect_size), rhs(vect_size); 
	// cout << "here ... 4" << endl;

	// the acceleration object
	converg_acceler<VectorXcd> Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update); 
	// cout << "here ... 5" << endl;
	// cout << "cond.maxStore = " << cond.maxStore << "; cond.wait = " << cond.wait << "; cond.rel_p = " << cond.rel_p << endl;
	// cout << "rhsBC =\n" << rhsBC << endl;
		   
 	// loop until f converges or until it's gone too long
	do { 
		// cout << "do: " << cts << endl;
		SC->bulkRHS_FE(cond, f, rhs, dummy);
		df = solver.solve(rhs+rhsBC) - f; // find the change in OP

		// if (method == string("acceleration"))
		Con_Acc.next_vector<MatrixXcd>( f, df, err ); // smart guess
		// else if (method == string("relaxation"))
		// f += 0.05*df;
		// err = df.norm()/f.norm();
		cts++;         // increment counter

		// output approx. percent completed
		cout << "\033[A\33[2K\r" << "\testimated: " << round((cts*100.)/cond.N_loop) << "% done; current error: " << err << endl;
	} while(err > cond.ACCUR && cts < cond.N_loop);

	if (err < cond.ACCUR) cout << "\tFound solution:" << endl;
	else cout << "\tResult did not converge satisfactorily:" << endl;
	cout << "\t\titerations = " << cts << endl;
	cout << "\t\trelative error = " << err << endl;

	return;
}
