//#include "includes.hpp"
#include "readwrite.hpp"

// most of the typedefs and custom structures are in this file:
#include "structures.hpp"
#include "SC_classes.hpp"

using namespace std;
using namespace Eigen;

// prototype for function defined in linear_eq_solver.cpp
void Solver(VectorXcd & f, SpMat_cd M, VectorXcd rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC);

int main()
{
	const int Nop=3; // number of the order parameter components

	// get all the information from the "conditions.txt"
	in_conditions cond;
	Bound_Cond eta_BC[Nop]; // boundary conditions for OP components
	Matrix2d **gradK;       // gradient coefficients in GL functional
	gradK = new Matrix2d *[Nop]; 
	for (int i = 0; i < Nop; i++) gradK[i] = new Matrix2d [Nop];

	read_input_data(Nop, cond, eta_BC, gradK, "conditions"+to_string(Nop)+".txt");
	// confirm_input_data(Nop, cond, eta_BC, gradK);
	
	// default parameters for the Convergence Accelerator
	cond.maxStore = 7; // 4
	cond.rel_p = 0.05; // 0.1
	cond.wait = 1;     // 2

	// if you want to change the values ... should we put these back into the conditions file?
	cout << "The default parameters are:\n\tmaxStore = " << cond.maxStore << "\n\trel_p = " << cond.rel_p << "\n\twait = " << cond.wait << endl;
	cout << "Do you want to change these values? (y/n): ";
	string res;
	cin >> res;
	if (res == "y" || res == "Y") {
		cout << "maxStore: ";
		cin >> cond.maxStore;
		cout << "rel_p: ";
		cin >> cond.rel_p;
		cout << "wait: ";
		cin >> cond.wait;
	}

	vector<int> no_update; // the vector of all the indeces of the OPvector that we don't want to change

	int GridSize = cond.SIZEu * cond.SIZEv; // number of grid points 
	int VectSize = cond.Nop * GridSize;     // number of elements in the overall OP vector 

	VectorXcd OPvector(VectSize);                // the vector of all values of the OP on the mesh
	VectorXcd rhsBC = VectorXcd::Zero(VectSize); // the addition to the rhs for D-type BC's
	VectorXcd dummy(VectSize);                   // a dummy vector, ...for free energy?
	SpMat_cd M(VectSize,VectSize);               // the matrix we will use to solve
	VectorXd freeEb(GridSize), freeEg(GridSize); // free energy on the grid points
	// add other observables here as desired...

	cout << "initializing object...";
	SC_class *pSC;
	// the SC object depending on given OP size
		 if (Nop == 3) pSC = new ThreeCompHe3( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
	else if (Nop == 5) pSC = new FiveCompHe3 ( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
	else {cout << "Unknown OP size. Exiting..." << endl; return 0;}
	cout << "done" << endl;

	cout << "initializing guess...";
	// initializeOPguess(cond, eta_BC, OPvector, GridSize, no_update); // set the OP vector to a good guess based on BC's
	pSC->initialOPguess(eta_BC, OPvector, no_update);
	cout << "done" << endl;

	// write the initial guess to file for debugging
	// WriteToFile(OPvector, "initGuess"+to_string(Nop)+".txt", Nop, cond);

	cout << "building solver matrix...";
	pSC->BuildSolverMatrix( M, rhsBC, OPvector, eta_BC, gradK );
	cout << "done" << endl;

	// cout << endl << "M =\n" << M << endl;

	cout << "solving system...";
	Solver(OPvector, M, rhsBC, cond, no_update, pSC);
	cout << "done" << endl;

	cout << "writing solution to file...";
	WriteToFile(OPvector, "solution"+to_string(Nop)+".txt", Nop, cond);
	cout << "done!" << endl;

	// ---- updated May 12, 2020 ----- 
	//* PLAN for the rest of the code:

	pSC->bulkRHS_FE(cond, OPvector, dummy, freeEb);
	// get the bulk contribution to free energy 
	WriteToFile(dummy, "bulkRHS_FE"+to_string(Nop)+".txt", Nop, cond);

	/* CONTINUE HERE!
	pSC->gradFE(freeEg, cond, OPvector, eta_BC, gradK);
	// get the gradient contribution to free energy 

	// then follows the print-out of the results into appropriate files
	// ... save the FE things as well!
	//*/

	//------  de-allocating gradK array ------------------
	for(int i = 0; i <Nop; i++) delete[] gradK[i]; 
	delete[] gradK; 
	delete pSC;

	return 0;
}