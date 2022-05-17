#include "SC_classes.hpp"


#include "structures.hpp"
#include <vector>

using namespace std;
using namespace Eigen;

// implement class specific
// CONTINUE HERE: TODO: implement this!!
// void gradFE(VectorXd & freeEg, in_conditions parameters, const T_vector OPvector,     Bound_Cond eta_BC[], Matrix2d **gradK);

// function to get the RHS of GL differential equations, determined from the bulk functional; 
// ------  coded in he3bulk.cpp ---------
// input: temp t, pressure p, matrix A, 
// output: betaB parameter that gives the bulk free energy of B-phase; derivative matrix dFdA, bulk free energy FE
void he3bulk(double t, double p, double & betaB, Eigen::Matrix<T_scalar,3,3> A, Eigen::Matrix<T_scalar,3,3> & dFdA, double & FE);

// --- Nop, grid_size are variables stored in base class

void ThreeCompHe3::bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, VectorXd & FEb) {
	// cout << "ThreeCompHe3::bulkRHS_FE" << endl;
	int grid_point;

	for (int v = 0; v < Nv; v++) 
	for (int u = 0; u < Nu; u++) {
		// cout << "v = " << v << endl;
		// cout << "u = " << u << endl;
		
		// we go over all grid points, 
		grid_point = ID(u,v,0);
		// cout << "grid_point = " << grid_point << endl;
		
		// and for all grid points collect all OP components into one eta_bulk[] array
		for ( int i = 0; i < Nop; i++) eta_bulk[i] = OPvector( ID(u, v, i) );
		// cout << "here 1" << endl;

		// create He3 OP matrix
		Matrix<T_scalar,3,3> A; 
					   A << eta_bulk[0],     0,         0,
						  0,          eta_bulk[1],    0, 
						  0,               0,     eta_bulk[2];
		// cout << "here 2" << endl;

		Matrix<T_scalar,3,3> dFdA;
		// cout << "here 3" << endl;
		double FEbulk, betaB;
		// cout << "here 4" << endl;
		T_scalar dFdeta[3];
		// cout << "here 5" << endl;

		// find the derivatives and bulk Free energy
		he3bulk(parameters.T, parameters.P, betaB, A, dFdA, FEbulk);
		// cout << "here 6" << endl;
		dFdeta[0] = dFdA(0,0);
		// cout << "here 7" << endl;
		dFdeta[1] = dFdA(1,1);
		// cout << "here 8" << endl;
		dFdeta[2] = dFdA(2,2);
		// cout << "here 9" << endl;

		// and put them into the big vector in the same order as OP components
		for ( int i = 0; i < Nop; i++) newRHSvector( ID(u, v, i) ) = dFdeta[i]; 
		// cout << "here 10" << endl;

		// cout << "FEbulk = " << FEbulk << endl;

		FEb(grid_point) = FEbulk;
		// cout << "here 11" << endl;
	}

	return;
}

void ThreeCompHe3::gradFE(Eigen::VectorXd & freeEg, const T_vector OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
	//// Follow equation 10 in Latex doc
	// we will modify the freeEg vector
	freeEg *= 0.; // make sure it's all 0
	VectorXcd temp;

	/*
	for (int n = 0; n < Nop; n++) {
		for (int m = 0; m < Nop; m++) {
			SpMat_cd lhs(grid_size,2);
			lhs(0,0) = this->Du_BD(eta_BC[m],m,temp,temp) * OPvector( seq(m*grid_size,(m+1)*grid_size) );
			lhs(0,1) = this->Dv_BD(eta_BC[m],m,temp,temp) * OPvector( seq(m*grid_size,(m+1)*grid_size) );

			SpMat_cd rhs(2,grid_size);
			rhs(0,0) = this->Du_BD(eta_BC[n],n,temp,temp) * OPvector( seq(n*grid_size,(n+1)*grid_size) );
			rhs(1,0) = this->Du_BD(eta_BC[n],n,temp,temp) * OPvector( seq(n*grid_size,(n+1)*grid_size) );

			freeEg += lhs * gradK[m][n] * rhs;
		}
	}
	*/

	// TODO: CONTINUE HERE:
	// // How can we put this all into a matrix?
	// for (int n = 0; n < Nop; n++) {
	// 	for (int m = 0; m < Nop; m++) {
			
	// 		auto opV_slice = OPvector( seq(m*grid_size, (m+1)*grid_size) );

	// 		Eigen::Matrix<T_scalar,-1,-1> lhs(grid_size,2);
	// 		lhs << Du_BD(eta_BC[m], m, temp, temp) * opV_slice,
	// 				Dv_BD(eta_BC[m], m, temp, temp) * opV_slice;

	// 		opV_slice = OPvector( seq(n*grid_size, (n+1)*grid_size) );

	// 		Eigen::Matrix<T_scalar,-1,-1> rhs(2,grid_size);
	// 		rhs << Du_BD(eta_BC[n], n, temp, temp) * opV_slice,
	// 				Dv_BD(eta_BC[n], n, temp, temp) * opV_slice;

	// 		freeEg += lhs.conjugate() * gradK[m][n] * rhs;
	// 	}
	// }
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
					  A <<  eta_bulk[0],    0,      eta_bulk[3],
					          0,          eta_bulk[1],    0,
						   eta_bulk[4],    0,      eta_bulk[2];

		Matrix<T_scalar,3,3> dFdA;
		double FEbulk, betaB;
		T_scalar dFdeta[5];
		// find the derivatives and bulk Free energy
		he3bulk(parameters.T, parameters.P, betaB, A, dFdA, FEbulk); 
		dFdeta[0] = dFdA(0,0);
		dFdeta[1] = dFdA(1,1);
		dFdeta[2] = dFdA(2,2);
		dFdeta[3] = dFdA(0,2);
		dFdeta[4] = dFdA(2,0);
		// and put them into the big vector in the same order as OP components 
		for ( int i = 0; i < Nop; i++) newRHSvector( ID(u, v, i) ) = dFdeta[i];
		FEb(grid_point) = FEbulk;
	}

	return;
}
void FiveCompHe3::gradFE(Eigen::VectorXd & freeEg, const T_vector OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
	//
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
void OneCompSC::gradFE(Eigen::VectorXd & freeEg, const T_vector OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
	//
}

