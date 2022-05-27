#include "SC_classes.hpp"


#include "structures.hpp"
#include <vector>

using namespace std;
using namespace Eigen;

// function to get the RHS of GL differential equations, determined from the bulk functional; 
// ------  coded in he3bulk.cpp ---------
// input: temp t, pressure p, matrix A, 
// output: betaB parameter that gives the bulk free energy of B-phase; derivative matrix dFdA, bulk free energy FE
void he3bulk(double t, double p, double & betaB, Eigen::Matrix<T_scalar,3,3> A, Eigen::Matrix<T_scalar,3,3> & dFdA, double & FE);

// --- Nop, grid_size are variables stored in base class

void ThreeCompHe3::bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, VectorXd & FEb) {
	int grid_point;

	for (int v = 0; v < Nv; v++) 
	for (int u = 0; u < Nu; u++) {
		// we go over all grid points, 
		grid_point = ID(u,v,0);
		
		// and for all grid points collect all OP components into one eta_bulk[] array
		for ( int i = 0; i < Nop; i++) eta_bulk[i] = OPvector( ID(u, v, i) );

		// create He3 OP matrix
		Matrix<T_scalar,3,3> A; 
					   A << eta_bulk[0],     0,         0,
						  0,          eta_bulk[1],    0, 
						  0,               0,     eta_bulk[2];

		Matrix<T_scalar,3,3> dFdA;
		double FEbulk, betaB;
		T_scalar dFdeta[3];

		// find the derivatives and bulk Free energy
		he3bulk(parameters.T, parameters.P, betaB, A, dFdA, FEbulk);
		dFdeta[0] = dFdA(0,0);
		dFdeta[1] = dFdA(1,1);
		dFdeta[2] = dFdA(2,2);

		// and put them into the big vector in the same order as OP components
		for ( int i = 0; i < Nop; i++) newRHSvector( ID(u, v, i) ) = dFdeta[i]; 

		FEb(grid_point) = FEbulk;
	}

	return;
}
void ThreeCompHe3::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
	// Follow equation 10 and 36 in Latex doc
	// cout << "In 'ThreeCompHe3::gradFE'" << endl;
	// freeEg = Eigen::VectorXd(grid_size*grid_size); // resize the free energy grad vector
	T_vector eta_dag = OPvector.adjoint();         // eta daggar; for more readable calculations
	T_vector F_times_eta = FEgrad * OPvector;      // the rhs vector in equation 36

	// loop over all points on the grid
	for (int u = 0; u < Nu; u++) {
		for (int v = 0; v < Nv; v++) {

			int FE_id = ID(u,v,0); // calculate the ID of the point
			freeEg( FE_id ) = 0;   // make sure the element starts at 0

			// loop over all OP comonents...get all contributions to the FE at this point
			for (int n = 0; n < Nop; n++) {
				int id = ID(u,v,n);
				freeEg( FE_id ) += ( eta_dag(id) * F_times_eta(id) ).real();
			}
		}
	}
	freeEg /= h*h;
	// cout << "freeEg = " << freeEg << endl;
}


void FiveCompHe3::bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, VectorXd & FEb) {
	int grid_point; 

	for (int v = 0; v < Nv; v++) 
	for (int u = 0; u < Nu; u++){ 
		// we go over all grid points, 
		grid_point = ID(u,v,0);
		// and for all grid points collect all OP components into one eta_bulk[] array
		for ( int i = 0; i < Nop; i++) eta_bulk[i] = OPvector( ID(u, v, i) ); 

		// create He3 OP matrix
		Matrix<T_scalar,3,3> A; 
					  A <<  eta_bulk[0],    0,      eta_bulk[4],
					          0,          eta_bulk[1],    0,
						   eta_bulk[3],    0,      eta_bulk[2];

		Matrix<T_scalar,3,3> dFdA;
		double FEbulk, betaB;
		T_scalar dFdeta[5];
		// find the derivatives and bulk Free energy
		he3bulk(parameters.T, parameters.P, betaB, A, dFdA, FEbulk); 
		dFdeta[0] = dFdA(0,0);
		dFdeta[1] = dFdA(1,1);
		dFdeta[2] = dFdA(2,2);
		dFdeta[3] = dFdA(2,0);
		dFdeta[4] = dFdA(0,2);
		// and put them into the big vector in the same order as OP components 
		for ( int i = 0; i < Nop; i++) newRHSvector( ID(u, v, i) ) = dFdeta[i];
		FEb(grid_point) = FEbulk;
	}

	return;
}
void FiveCompHe3::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
	// freeEg = Eigen::VectorXd(grid_size*grid_size);
	T_vector eta_dag = OPvector.adjoint();
	T_vector F_times_eta = FEgrad * OPvector;
	for (int u = 0; u < Nu; u++) {
		for (int v = 0; v < Nv; v++) {
			int FE_id = ID(u,v,0);
			freeEg( FE_id ) = 0;
			for (int n = 0; n < Nop; n++) {
				int id = ID(u,v,n);
				freeEg( FE_id ) += ( eta_dag(id) * F_times_eta(id) ).real();
			}
		}
	}
	freeEg /= h*h;
}
double FiveCompHe3 :: defectEnergy(const Eigen::VectorXd & freeEb, const Eigen::VectorXd & freeEg) {
	double DefectEnergy = 0.;
	for (int u = 0; u < Nu; u++) {
		for (int v = 0; v < Nv; v++) {
			int FE_id = ID(u,v,0);
			DefectEnergy += freeEg(FE_id);
		}
	}
	return DefectEnergy/(h*h);
}


void OneCompSC::bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, VectorXd & FEb) {
	for (int grid_point = 0; grid_point < grid_size; grid_point++){ 
		eta_bulk[0] = OPvector( grid_point ); 

		double FEbulk;
		T_scalar dFdeta;

		// FE_bulk = -alpha*eta^2 + beta eta^4
		// find the derivatives and bulk Free energy
		dFdeta = -eta_bulk[0] +norm(eta_bulk[0])*eta_bulk[0];
		FEbulk = -2.0*norm(eta_bulk[0]) + pow(norm(eta_bulk[0]),2); 
		newRHSvector( grid_point ) = dFdeta;
		FEb(grid_point) = FEbulk;
	}

	return;
}
void OneCompSC::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
	//
}

