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
	// for(int i=0; i< rhs_local.size(); i++) if(abs(rhs_local[i])==0) rhs_local[i]=h2_times_rhs_bulk[i];
	rhs_local = h2_times_rhs_bulk + rhsBC;
	return rhs_local;
}

void Solver(T_vector & f, SpMat_cd M, T_vector rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC, string method)
{
	int grid_size=cond.SIZEu * cond.SIZEv; 
	int vect_size=cond.Nop * grid_size; 
	double h2 = cond.STEP*cond.STEP;
	// FILE *out; out=fopen("iterations.dat","w");

	if (no_update.size()) cout << "using no_update" << endl;

	// use LU decomposition to solve the system
	SparseLU< SpMat_cd, COLAMDOrdering<int> > solver;
   	solver.analyzePattern( M ); 
	solver.factorize( M );

	// check to see if the solver failed
	if (solver.info() == Eigen::Success) cout << "\tSolver: successfully built" << endl;
	else if (solver.info() == Eigen::NumericalIssue) { // for debugging non-invertable matrices
		cout << "Solver: numerical issues" << endl;
		// for (int k=0; k < SolverMatrix.outerSize(); ++k)
		// 	for (SparseMatrix<Scalar_type>::InnerIterator it(SolverMatrix,k); it; ++it)
		// 		cout << "(" << it.row() << "," << it.col() << ")\t";
		return;
	}

	int cts = 0; // count loops
	double err;  // to store current error
	VectorXd dummy(vect_size); // dummy free energy variable for RHS function
	T_vector df(vect_size), rhs(vect_size); // 
	converg_acceler<T_vector> Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update); // the acceleration object
		   
 	// loop until f converges or until it's gone too long
	do { 
		/*
		for(int i=0; i<cond.SIZEv; i++) 
			fprintf(out, "%f  %f  %f  %f\n", i*cond.STEP, real(f[i+cond.SIZEv*0]), real(f[i+cond.SIZEv*1]), real(f[i+cond.SIZEv*2]));
		fprintf(out, "\n"); 
		*/

		//theoretical3comp_rhs(cond, f, rhs_th, dummy);
		// cout << "do: " << cts << endl;


		// save output of each guess for debugging
		// if (debug) {
		// 	// ...
		// 	WriteToFile(f, "debugging_files/solu_"+to_string(cond.Nop)+"iter_"+to_string(cts)+".txt", cond);
		// }

		SC->bulkRHS_FE(cond, f, rhs, dummy);
		df = solver.solve( get_rhs(h2*rhs, rhsBC) ) - f; // find the change in OP

		// acceleration method
		if (method == "acceleration") Con_Acc.next_vector<T_matrix>( f, df, err ); // smart guess
		else { // normal relaxation method
			f += cond.rel_p*df;
			err = df.norm()/f.norm();
		}

		cts++;         // increment counter
		// output approx. percent completed
		//cout << endl;
		//cout << "\033[A\33[2K\r" << "\testimated: " << round((cts*100.)/cond.N_loop) << "% done; current error: " << err << endl;
	} while(err > cond.ACCUR && cts < cond.N_loop);

	if (err < cond.ACCUR) cout << "\tFound solution:" << endl;
	else cout << "\tResult did not converge satisfactorily:" << endl;
	cout << "\t\titerations = " << cts << endl;
	cout << "\t\trelative error = " << err << endl;

	return;
}

/*
  // --------------- some test functions ------------------------------------------
void theoretical3comp_rhs(in_conditions cond, T_vector & OPvector, T_vector & RHSvector, Eigen::VectorXd & FEb);

void theoretical3comp_rhs(in_conditions cond, T_vector & OPvector, T_vector & RHSvector, Eigen::VectorXd & FEb)
{
	complex<double> Axx, Ayy, Azz;
	complex<double> Rxx, Ryy, Rzz;

	for(int i=0; i<cond.SIZEv; i++){
			Axx= OPvector[i+cond.SIZEv*0]; 
			Ayy= OPvector[i+cond.SIZEv*1]; 
			Azz= OPvector[i+cond.SIZEv*2]; 
			Rxx= -Axx + pow(Axx,3) - 1.0/5.0 * (Axx*Axx - Azz*Azz)*Axx;
			Ryy=Rxx;
			Rzz= -Azz + pow(Azz,3) + 2.0/5.0 * (Axx*Axx - Azz*Azz)*Azz;
			RHSvector[i+cond.SIZEv*0]=Rxx;
			RHSvector[i+cond.SIZEv*1]=Ryy;
			RHSvector[i+cond.SIZEv*2]=Rzz;
	}

	return;
}
*/
