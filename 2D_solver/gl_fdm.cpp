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

	// get all the information from the conditions<...>.txt
	in_conditions cond;
	vector<Bound_Cond> eta_BC; // boundary conditions for OP components

	read_input_data(Nop, cond, eta_BC, file_name);
	// confirm_input_data(Nop, cond, eta_BC); // confirm the input by printing it out

	vector<int> no_update; // the vector of all the indeces of the OPvector that we don't want to change

	int GridSize = cond.SIZEu * cond.SIZEv; // number of grid points 
	int VectSize = cond.Nop * GridSize;     // number of elements in the overall OP vector

	T_vector OPvector(VectSize);               // the vector of all values of the OP on the mesh
	T_vector rhsBC = T_vector::Zero(VectSize); // the addition to the rhs for D-type BC's
	T_vector dummy(VectSize);                  // a dummy vector, ...for free energy?
	SpMat_cd M(VectSize,VectSize);             // the matrix we will use to solve
	VectorXd FEdens(GridSize), freeEb(GridSize), freeEg(GridSize); // free energy on the grid points
	// add other observables here as desired...

	// ===============================================================================================================

	cout << "initializing object and guess..." << endl;
	SC_class *pSC; // the SC object...
	// ... depending on given OP size
	if (argc == 3 && *(argv[2]) == 'c') { // if it is for a cylindrical system
		pSC = new Cylindrical ( Nop, cond.SIZEu, cond.SIZEv, cond.STEP, eta_BC );
		if (Nop == 5) {
			if (file_name == string("conditions5c.txt"))              pSC->initialOPguess_Cylindrical_simple5(eta_BC, OPvector, no_update);
			else if (file_name == string("conditions5c_AzzFlip.txt")) pSC->initialOPguess_Cylindrical_AzzFlip(eta_BC, OPvector, no_update);
			else if (file_name == string("conditions5c_bubble.txt"))  pSC->initialOPguess_Cylindrical_bubble (eta_BC, OPvector, no_update);
			else {
				cout << "Unknown file_name. Exiting..." << endl;
				delete pSC;
				return 1;
			}
		} else if (Nop == 3) {
			if (file_name == string("conditions3c.txt")) pSC->initialOPguess_Cylindrical_simple3(eta_BC, OPvector, no_update);
			else {
				cout << "WARNING: you either added the optional argument '[c]' on accident, or passed in the wrong file name for a cylindrical system." << endl;
				delete pSC;
				return 1;
			}
		}
	}
	// other wise we'll use the cartesian system
	else if (Nop == 3) {
		pSC = new ThreeCompHe3( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
		if (file_name == string("conditions3.txt")) pSC->initialOPguess(eta_BC, OPvector, no_update); // set the OP vector to a good guess based on BC's
		else {
			cout << "WARNING: you probably forgot the optional argument '[c]'." << endl;
			delete pSC;
			return 1;
		}
	} else if (Nop == 5) {
		pSC = new FiveCompHe3( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
		// ---- set the OP vector to a good guess based on BC's ----
		if (file_name == string("conditions5.txt"))              pSC->initialOPguess             (eta_BC, OPvector, no_update);
		else if (file_name == string("conditions5_AyyFlip.txt")) pSC->initGuessWithCircularDomain(eta_BC, OPvector, no_update);
		else if (file_name == string("conditions5_wall.txt"))    pSC->initOPguess_special(cond, eta_BC, OPvector, no_update); // Anton's version
		else {
			cout << "Unknown file_name. Exiting..." << endl;
			delete pSC;
			return 1;
		}
	} else if (Nop == 1) {
		pSC = new OneCompSC( Nop, cond.SIZEu, cond.SIZEv, cond.STEP );
	} else {
		cout << "Unknown OP size. Exiting..." << endl;
		delete pSC; // IS THIS ACTUALLY OK IF pSC HAS NOT BEEN SET?
		return 1;
	}
	
	// ===============================================================================================================

	cout << "building solver matrix...";
	if (argc == 3 && *(argv[2]) == 'c') // TODO: combine these functions into one?
		pSC->BuildSolverMatrixCyl( M, rhsBC, OPvector, eta_BC );
	else
		pSC->BuildSolverMatrix( M, rhsBC, OPvector, eta_BC );
	cout << "done" << endl;

	// ===============================================================================================================

	cout << "solving system..." << endl;
	Solver(OPvector, M, rhsBC, cond, no_update, pSC); // solve the system setup above
	cout << "solved!" << endl;

	cout << "calculating Free energy ...";
	double totalFE = pSC->FreeEn(OPvector, cond, FEdens, freeEb, freeEg);
	cout << "done" << endl;
	cout << "the total energy loss in DW is " << totalFE << "\n";

	// ===============================================================================================================

	// write everything to file
	pSC->WriteAllToFile(OPvector, freeEb, freeEg, "output_OP"+to_string(Nop)+( (argc == 3 && *(argv[2]) == 'c') ? "c" : "" )+".txt");

	// ===============================================================================================================

	// de-allocating memory
	delete pSC;

	return 0;
}
