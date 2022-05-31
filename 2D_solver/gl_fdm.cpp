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
	bool debug = false;
	if (argc == 3) {
		// using debugging!
		cout << "NOTICE: using debugging methods...\nDo you want to proceed? (y/n) ";
		string ans;
		cin >> ans;
		if (ans=="y"||ans=="Y") debug = true;
		else cout << endl << "NOT USING DEBUGGING ANYMORE." << endl;
	}
	else if (argc != 2) {cout << "ERROR: need an argument for 'Nop'; do so like: '$ ./gl_fdm 3'." << endl; return 1;}
	const int Nop = *(argv[1]) - '0'; // read in the int from the terminal call

	// the list of K-matrix components to be used in building the small K-matrix
	vector<bool> components;
	components.insert(components.end(),{true,true,true,true});
	// build the small K-matrix based on the OP size
	auto smKMat = kMatrix::smallKMatrix(Nop, components);

	// get all the information from the "conditions.txt"
	in_conditions cond;
	Bound_Cond eta_BC[Nop];      // boundary conditions for OP components

	read_input_data(Nop, cond, eta_BC, "conditions"+to_string(Nop)+".txt");
	if (debug) confirm_input_data(Nop, cond, eta_BC, smKMat);
	
	// default parameters for the Convergence Accelerator
	cond.maxStore = 5; // 4
	cond.rel_p = 0.1;  // 0.1
	cond.wait = 1;     // 2

	// if you want to change the values ... should we put these back into the conditions file?
	if (debug) {
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
	}

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

	cout << "initializing object...";
	SC_class *pSC; // the SC object...
	// ... depending on given OP size
		 if (Nop == 3) pSC = new ThreeCompHe3( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
	else if (Nop == 5) pSC = new FiveCompHe3 ( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
	else if (Nop == 1) pSC = new OneCompSC   ( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
	else {cout << "Unknown OP size. Exiting..." << endl; return 0;}
	cout << "done" << endl;

	// ===============================================================================================================

	cout << "initializing guess...";
	// - - - - switch these here to initialize in other ways
	pSC->initialOPguess(eta_BC, OPvector, no_update); // set the OP vector to a good guess based on BC's
	// pSC->initOPguess_special(cond, eta_BC, OPvector, gradK, no_update); // Anton's version
	// pSC->initGuessWithCircularDomain(eta_BC, OPvector, no_update); // CONTINUE HERE! TESTING THIS ONE!
	cout << "done" << endl;

	if (debug) { // write the initial guess to file, for debugging
		pSC->WriteToFile(OPvector, "initGuess"+to_string(Nop)+".txt", 1);
	}

	// ===============================================================================================================

	cout << "building solver matrix...";
	pSC->BuildSolverMatrix( M, rhsBC, OPvector, eta_BC, smKMat );
	cout << "done" << endl;

	// if (debug) { // For debugging only...shouldn't print if gsize > ~10^2
	// 	cout << endl << "M =\n" << M << endl;
	// }

	// ===============================================================================================================

	cout << "solving system..." << endl;
	Solver(OPvector, M, rhsBC, cond, no_update, pSC); // solve the system setup above
	cout << "solved!" << endl;

	cout << "calculating bulkRHS_FE...";
	pSC->bulkRHS_FE(cond, OPvector, dummy, freeEb); // get the bulk contribution to free energy
	cout << "done" << endl;

	cout << "calculating gradFE...";
	pSC->gradFE(freeEg, OPvector, eta_BC, smKMat); // get the gradient contribution to free energy
	cout << "done" << endl;

	// ===============================================================================================================

	// write everything to file
	pSC->WriteAllToFile(OPvector, freeEb, freeEg, "output_OP"+to_string(Nop)+".txt");

	// calculate the defect energy
	// cout << "Energy defect: " << pSC->defectEnergy(freeEb, freeEg) << endl;

	// ===============================================================================================================

	//------  de-allocating smKMat array ----------
	for(int i = 0; i <Nop; i++) delete[] smKMat[i]; 
	delete[] smKMat; 
	delete pSC;

	return 0;
}
