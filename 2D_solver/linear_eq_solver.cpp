#include <iostream>
#include <fstream> // for file in/out
#include <string>
#include <complex>
#include <chrono> // for timing
#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>

#include "ConvergenceAccelerator.hpp"
#include "structures.hpp"
#include "readwrite.hpp"
#include "SC_classes.hpp"

using namespace Eigen;
using namespace std;

// solves self-consistency for linear system of equations
//     M*f = RHS(f)
// RHS function is given through *SC class
// on input f is initial guess, 
// on output it is the required solution 
//---------------------------------------------------------------------------------

T_vector get_rhs(T_vector h2_times_rhs_bulk, T_vector rhsBC){
	auto rhs_local=rhsBC;
	rhs_local = h2_times_rhs_bulk + rhsBC;
	return rhs_local;
}

void Solver(T_vector & f, const SpMat_cd & M, T_vector rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC) {
	int grid_size=cond.SIZEu * cond.SIZEv;
	int vect_size=cond.Nop * grid_size;
	double h2 = cond.STEP*cond.STEP;

	if (no_update.size()) cout << "using no_update" << endl;

	// use LU decomposition to solve the system
	SparseLU< SpMat_cd, COLAMDOrdering<int> > solver;
   	solver.analyzePattern( M ); 
	solver.factorize( M );

	// check to see if the solver failed
	if (solver.info() == Eigen::Success) cout << "\tSolver: successfully built" << endl;
	else if (solver.info() == Eigen::NumericalIssue) { // for debugging non-invertable matrices
		cout << "Solver: numerical issues" << endl;
		return;
	}

	int iter = 0; // count iterations 
	double err;   // to store current error
	VectorXd dummy(vect_size); // dummy free energy variable for RHS function
	T_vector df(vect_size), rhs(vect_size); 

	// the acceleration object
	converg_acceler<T_vector> Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update);
		   
 	// loop until f converges or until it's gone too long
	do {
		iter++; // increment counter

		SC->bulkRHS_FE(cond, f, rhs, dummy);
		df = (M*f - get_rhs(h2*rhs, rhsBC) ); // the best, works with both relaxation and acceleration methods 
		for(auto id : no_update) df(id) = 0.0; // no update for Dirichlet bc and other fixed points
		Con_Acc.next_vector<T_matrix>( f, df, err ); // smart next guess

		// output approximate amount completed
		if (iter%25==0) cout << "\033[A\33[2K\r" << "\t done " << iter << " iterations (out of " << cond.N_loop << ")" << "; current error: " << err << endl;
	} while(err > cond.ACCUR && iter < cond.N_loop);

	if (err < cond.ACCUR) cout << "\tFound solution:" << endl;
	else cout << "\tResult did not converge satisfactorily:" << endl;
	cout << "\t\titerations = " << iter << endl;
	cout << "\t\trelative error = " << err << endl;

	cout << endl;
	return;
}
