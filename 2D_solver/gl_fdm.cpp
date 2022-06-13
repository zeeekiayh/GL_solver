#include "readwrite.hpp"

// most of the typedefs and custom structures are in this file:
#include "structures.hpp"
#include "SC_classes.hpp"

using namespace std;
using namespace Eigen;

// ===========================================================
// Function prototypes defined in other files
	// prototype for function defined in linear_eq_solver.cpp
	void Solver(T_vector & f, SpMat_cd M, T_vector rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC);
// ===========================================================

int main(int argc, char** argv)
{
	int Nop;
	if (argc < 2 || argc > 3) {
		cout << "ERROR: need an argument for 'file_name'; do so like: '$ ./gl_fdm <file_name> [c]'." << endl;
		return 1;
	}
	string file_name = string(argv[1]);

	// get all the information from the "conditions.txt"
	in_conditions cond;
	Bound_Cond *eta_BC; // boundary conditions for OP components

	read_input_data(Nop, cond, eta_BC, file_name);
	// confirm_input_data(Nop, cond, eta_BC); // DO WE WANT TO PASS ", gK" AS AN ARGUMNET, EVER?

	vector<int> no_update; // the vector of all the indeces of the OPvector that we don't want to change

	int GridSize = cond.SIZEu * cond.SIZEv; // number of grid points 
	int VectSize = cond.Nop * GridSize;     // number of elements in the overall OP vector

	T_vector OPvector(VectSize);                 // the vector of all values of the OP on the mesh
	T_vector rhsBC = T_vector::Zero(VectSize);   // the addition to the rhs for D-type BC's
	T_vector dummy(VectSize);                    // a dummy vector, ...for free energy?
	SpMat_cd M(VectSize,VectSize);               // the matrix we will use to solve
	VectorXd freeEb(GridSize), freeEg(GridSize); // free energy on the grid points
	// add other observables here as desired...

	// ===============================================================================================================

	cout << "initializing object and guess..." << endl;
	SC_class *pSC; // the SC object...
	// ... depending on given OP size
	if (argc == 3 && *(argv[2]) == 'c') { // if it is for a cylindrical system
		pSC = new Cylindrical ( Nop, cond.SIZEu, cond.SIZEv, cond.STEP, eta_BC );
		if (Nop == 5) {
			if (file_name == string("conditions5c.txt"))         pSC->initialOPguess_Cylindrical_simple5(eta_BC, OPvector, no_update);
			if (file_name == string("conditions5c_AzzFlip.txt")) pSC->initialOPguess_Cylindrical_AzzFlip(eta_BC, OPvector, no_update); // TODO: make file, "conditions5c_AzzFlip.txt"
			if (file_name == string("conditions5c_bubble.txt"))  pSC->initialOPguess_Cylindrical_bubble (eta_BC, OPvector, no_update);
		} else if (Nop == 3) {
			// code this next!
			if (file_name == string("conditions3c.txt")) pSC->initialOPguess_Cylindrical_simple3(eta_BC, OPvector, no_update);
			else {cout << "WARNING: you either added the optional argument '[c]' on accident, or passed in the wrong file name for a cylindrical system." << endl; return 1;}
		}
	}
	// other wise we'll use the cartesian system
	else if (Nop == 3) {
		pSC = new ThreeCompHe3( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
		if (file_name == string("conditions3.txt")) pSC->initialOPguess(eta_BC, OPvector, no_update); // set the OP vector to a good guess based on BC's
		else {cout << "WARNING: you probably forgot the optional argument '[c]'." << endl; return 1;}
	} else if (Nop == 5) {
		pSC = new FiveCompHe3( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
		// ---- set the OP vector to a good guess based on BC's ----
		if (file_name == string("conditions5.txt"))         pSC->initialOPguess             (eta_BC, OPvector, no_update);
		if (file_name == string("conditions5_AyyFlip.txt")) pSC->initGuessWithCircularDomain(eta_BC, OPvector, no_update);
		if (file_name == string("conditions5_wall.txt"))    pSC->initOPguess_special(cond, eta_BC, OPvector, no_update); // Anton's version
	} else if (Nop == 1) {
		pSC = new OneCompSC( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
	} else {
		cout << "Unknown OP size. Exiting..." << endl;
		delete pSC;
		return 1;
	}

	// ===============================================================================================================

	// if (debug) { // write the initial guess to file, for debugging
	// 	// pSC->WriteToFile(OPvector, "initGuess"+to_string(Nop)+".txt", 1);
	// 	pSC->WriteAllToFile(OPvector, freeEb, freeEg, "initGuess_output_OP"+to_string(Nop)+".txt");
	// }

	// ===============================================================================================================

	cout << "building solver matrix...";
	if (argc == 3 && *(argv[2]) == 'c')
		pSC->BuildSolverMatrixCyl( M, rhsBC, OPvector, eta_BC );
	else
		pSC->BuildSolverMatrix( M, rhsBC, OPvector, eta_BC );
	cout << "done" << endl;

	// if (debug) { // For debugging only...shouldn't print if gsize > ~10^2
	// 	if (VectSize < 200)
	// 		cout << endl << "M =\n" << M << endl;
	// 	else
	// 		cout << "VectSize = " << VectSize << endl;
	// }
	// SparseMatrix<double> M(rows,cols);
	
	// cout << "Matrix:" << endl;
	// for (int k=0; k<M.outerSize(); ++k)
	// for (SpMat_cd::InnerIterator it(M,k); it; ++it) {
	// 	if (real(it.value()) > 5e3) {
	// 		cout << endl << "\tvalue: " << it.value() << endl;
	// 		cout << "\tposition: (" << it.row() << "," << it.col() << ")" << endl;
	// 	}
	// }
	// cout << "==============================" << endl << endl;

	// cout << "OPvector:" << endl;
	// cout << OPvector << endl;
	// // for (int i=0; i<OPvector.size(); i++)
	// // 	if (real(OPvector(i)) > 5e3)
	// // 		cout << "\tOPvector(" << i << ") = " << OPvector(i) << endl << endl;
	// cout << "==============================" << endl << endl;

	// cout << "rhsBC:" << endl;
	// cout << rhsBC << endl;
	// // for (int i=0; i<rhsBC.size(); i++)
	// // 	if (real(rhsBC(i)) > 5e3)
	// // 		cout << "\tOPvector(" << i << ") = " << rhsBC(i) << endl << endl;
	// cout << "==============================" << endl << endl;

	// ===============================================================================================================

	cout << "solving system..." << endl;
	Solver(OPvector, M, rhsBC, cond, no_update, pSC); // solve the system setup above
	cout << "solved!" << endl;

	cout << "calculating bulkRHS_FE...";
	pSC->bulkRHS_FE(cond, OPvector, dummy, freeEb); // get the bulk contribution to free energy
	cout << "done" << endl;

	cout << "calculating gradFE...";
	pSC->gradFE(freeEg, OPvector, eta_BC); // MAY NEED TO PASS IN ", gK" AS ANOTHER ARGUMENT? // get the gradient contribution to free energy
	cout << "done" << endl;

	// ===============================================================================================================

	// write everything to file
	pSC->WriteAllToFile(OPvector, freeEb, freeEg, "output_OP"+to_string(Nop)+".txt");

	// calculate the defect energy
	// cout << "Energy defect: " << pSC->defectEnergy(freeEb, freeEg) << endl;

	// ===============================================================================================================

	// de-allocating memory
	delete pSC;

	return 0;
}
