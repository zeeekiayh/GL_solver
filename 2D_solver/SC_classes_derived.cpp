#include "SC_classes.hpp"
#include "structures.hpp"
#include <vector>
#include <eigen/unsupported/Eigen/KroneckerProduct> // for the Kronecker product in Place_subMatrix

using namespace std;
using namespace Eigen;

// function to get the RHS of GL differential equations, determined from the bulk functional; 
// ------  coded in he3bulk.cpp ---------
// input: temp t, pressure p, matrix A, 
// output: betaB parameter that gives the bulk free energy of B-phase; derivative matrix dFdA, bulk free energy FE
void he3bulk(double t, double p, double & betaB, Eigen::Matrix<T_scalar,3,3> A, Eigen::Matrix<T_scalar,3,3> & dFdA, double & FE);

// prototype for function defined in 'SC_classes_base.cpp'
SpMat_cd Place_subMatrix(int i, int j, int size, SpMat_cd sm);

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
		dFdeta[3] = dFdA(0,2);
		dFdeta[4] = dFdA(2,0);
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

void Cylindrical::Build_curvilinear_matrices(Bound_Cond eta_BC[]) {
	I.resize(vect_size,vect_size); // set the identity matrix to the right size
	I.setIdentity();
	
	Dr.resize(vect_size,vect_size);
	Dz.resize(vect_size,vect_size);
	Dphi.resize(vect_size,vect_size);
	r_inv.resize(grid_size,grid_size);

	T_vector initOPvector; // ? should it actually be passed in as an argument in this function?
	T_vector rhsBC; // ? should it actually be passed in as an argument?
	
	// SpMat_cd R(grid_size,grid_size); // there is no need for the vect_size x vect_size r_inv?
	for (int u = 0; u < Nu; u++) {
		for (int v = 0; v < Nv; v++) {
			int id = ID(u,v,0);
			r_inv.coeffRef(id,id) = (u == 0) ? 1e10 : 1./(u*h);
		}
	}

	for (int n = 0; n < Nop; n++) {
		Dr += Place_subMatrix(n,n,vect_size,Du_BD(eta_BC[n],n,initOPvector,rhsBC));
		Dz += Place_subMatrix(n,n,vect_size,Dv_BD(eta_BC[n],n,initOPvector,rhsBC));
		// r_inv += Place_subMatrix(n,n,vect_size,R);
	}

	MatrixXd temp(9,9);
	temp << 0, 0, 0, 0, 0,-1,-1, 0, 0,
			0, 0, 0, 0, 0, 1, 1, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0,-1,
			0, 0, 0, 0, 0, 0, 0,-1, 0,
			1,-1, 0, 0, 0, 0, 0, 0, 0,
			1,-1, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 1, 0, 0, 0, 0,
			0, 0, 0, 1, 0, 0, 0, 0, 0;
	SpMat_cd small_I(grid_size,grid_size);
	small_I.setIdentity();
	Dphi = kroneckerProduct(temp,small_I);
}
void Cylindrical::BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
	// make all the matrices in eq. (54)
	SpMat_cd Dr2t(grid_size,grid_size); // D_r^2 ~
	SpMat_cd Drp(grid_size,grid_size);  // D_r^+
	SpMat_cd Drm(grid_size,grid_size);  // D_r^-
	SpMat_cd Drzp(grid_size,grid_size); // D_rz^+
	SpMat_cd D2t(grid_size,grid_size);  // D^2 ~

	// define D matrices that are used in others' definitions
	SpMat_cd Dr2_gsize(grid_size,grid_size),
			 Dr_gsize(grid_size,grid_size),
			 Dr_over_r(grid_size,grid_size),
			 Drz_gsize(grid_size,grid_size),
			 Dz_over_r(grid_size,grid_size),
			 Dz2_gsize(grid_size,grid_size),
			 Dz_gsize(grid_size,grid_size);
	SpMat_cd r2_inv = r_inv*r_inv; // this one is the same regardless of BC's

	// the matrices -- 2 parts of the solver matrix
	SpMat_cd DK23(vect_size,vect_size), DK1(vect_size,vect_size);

	for (int n = 0; n < Nop; n++) {
		// initialize the D matrices based on BC's
		// only make those that we ABOLUTELY need, because these take a LONG time to build!
		Dr2_gsize = Du2_BD(eta_BC[n],n,initOPvector,rhsBC);
		Dz2_gsize = Dv2_BD(eta_BC[n],n,initOPvector,rhsBC);
		if (n == 0 || n == 2 || n == 3 || n == 4) Drz_gsize = Duv_BD(eta_BC[n],n,initOPvector,rhsBC);
		if (n == 0 || n == 1) Dr_gsize  = Du_BD (eta_BC[n],n,initOPvector,rhsBC);
		if (n == 0) Dz_gsize  = Du_BD (eta_BC[n],n,initOPvector,rhsBC);
		Dr_over_r = Dr_gsize*r_inv;
		if (n == 0 || n == 1 || n == 3 || n == 4) Dz_over_r = Dz_gsize*r_inv;

		// build each complicated D matrix
		if (n == 0 || n == 3) Dr2t = Dr2_gsize + Dr_over_r - r2_inv;
		if (n == 0) Drp = Dr_gsize + r_inv;
		if (n == 1) Drm = Dr_gsize - r_inv;
		if (n == 0 || n == 3) Drzp = Drz_gsize + Dz_over_r;
		D2t = Dz2_gsize + Dr2_gsize + Dr_over_r - r2_inv;

		// TODO/WARNING: this is hard-coded for the 5-component system!!
		// Can we, or do we even need to make it for larger systems as well?

		// build eq.'s (55-56)
		if (n == 0) {
			DK23 += Place_subMatrix(0,n,vect_size,Dr2t)
				  + Place_subMatrix(1,n,vect_size,Drp*r_inv)
				  + Place_subMatrix(4,n,vect_size,Drzp);
			
			DK1 += Place_subMatrix(0,n,vect_size,D2t-r2_inv)
			     + Place_subMatrix(1,n,vect_size,2.*r2_inv);
		} else if (n == 1) {
			DK23 += Place_subMatrix(0,n,vect_size,-Drm)
				  + Place_subMatrix(1,n,vect_size,-r2_inv)
				  + Place_subMatrix(4,n,vect_size,-Dz_over_r);
			
			DK1 += Place_subMatrix(1,n,vect_size,D2t-r2_inv)
			     + Place_subMatrix(0,n,vect_size,2.*r2_inv);
		} else if (n == 2) {
			DK23 += Place_subMatrix(2,n,vect_size,Dz2_gsize)
			      + Place_subMatrix(3,n,vect_size,Drz_gsize);
			
			DK1 += Place_subMatrix(n,n,vect_size,D2t+r2_inv);
		} else if (n == 3) {
			DK23 += Place_subMatrix(2,n,vect_size,Drzp)
			      + Place_subMatrix(3,n,vect_size,Dr2t);
			
			DK1 += Place_subMatrix(n,n,vect_size,D2t);
		} else if (n == 4) {
			DK23 += Place_subMatrix(0,n,vect_size,Drz_gsize)
				  + Place_subMatrix(1,n,vect_size,Dz_over_r)
				  + Place_subMatrix(4,n,vect_size,Dz2_gsize);
			
			DK1 += Place_subMatrix(n,n,vect_size,D2t);
		}
	}

	// constants
	int K1 = 1, K23 = 2;
	M = K1 * DK1 + K23 * DK23;
}
void Cylindrical::bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & freeEb) {
	//
}
void Cylindrical::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
	//
}
double Cylindrical::defectEnergy(const Eigen::VectorXd & freeEb, const Eigen::VectorXd & freeEg) {
	//
}
