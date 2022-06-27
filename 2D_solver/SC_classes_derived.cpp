#include "SC_classes.hpp"
#include "readwrite.hpp"
#include "structures.hpp"
#include <eigen/unsupported/Eigen/KroneckerProduct> // for the Kronecker product in Build_curvilinear_matrices

using namespace std;
using namespace Eigen;

// ===========================================================
// Function prototypes defined in other files
	// function to get the RHS of GL differential equations, determined from the bulk functional; 
	// ------  coded in he3bulk.cpp ---------
	// input: temp t, pressure p, matrix A, 
	// output: betaB parameter that gives the bulk free energy of B-phase; derivative matrix dFdA, bulk free energy FE
	void he3bulk(double t, double p, double & betaB, Eigen::Matrix<T_scalar,3,3> A, Eigen::Matrix<T_scalar,3,3> & dFdA, double & FE);

	// For building the solver matrices
	// ------  coded in SC_classes_base.cpp ---------
	SpMat_cd Place_subMatrix(int i, int j, int size, SpMat_cd sm);

	// For solving some systems as an inital guess for others
	// ------  coded in linear_eq_solver.cpp ---------
	void Solver(T_vector & f, SpMat_cd M, T_vector rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC);
// ===========================================================


// ===========================================================
// ThreeCompHe3 :: function definitions
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

	// on input FEdens is the FE density of a reference configuration 
	double ThreeCompHe3::FreeEn(T_vector & OPvector, in_conditions cond, Eigen::VectorXd & FEdens, Eigen::VectorXd & FEb, Eigen::VectorXd & FEg) {
		T_vector dummy(vect_size); 
		this->bulkRHS_FE(cond, OPvector, dummy, FEb); // get the bulk contribution to free energy density

		T_vector Du_eta(vect_size), Dv_eta(vect_size); 
		if( Nu > 1 ) Du_eta = Du_FE * OPvector; else Du_eta = T_vector::Zero(vect_size);
		if( Nv > 1 ) Dv_eta = Dv_FE * OPvector; else Dv_eta = T_vector::Zero(vect_size);

		T_vector Du_eta_dag=Du_eta.adjoint(), Dv_eta_dag=Dv_eta.adjoint(); 
		
		double FE_minus_ref=0.0; // relative to some reference free energy 
		int x = 0, z = 1;            // indices for the K-matrix
		double wx, wz;               // weights for the integral free energy (trapezoidal rule)

		for (int u = 0; u < Nu; u++) {
			if( (u==0 || u==Nu-1) && Nu>1 ) wx=0.5; else wx=1.0;
			for (int v = 0; v < Nv; v++) {
				if( (v==0 || v==Nv-1) && Nv>1 ) wz=0.5; else wz=1.0;

				int id = ID(u,v,0);
				FEg( id ) = 0;

				for (int m = 0; m < Nop; m++) {
					for (int n = 0; n < Nop; n++) {
						int id_m = ID(u,v,m);
						int id_n = ID(u,v,n);

						FEg( id ) += gK[m][n](x,x) * real(Du_eta_dag(id_m) * Du_eta(id_n)); 
						FEg( id ) += gK[m][n](z,z) * real(Dv_eta_dag(id_m) * Dv_eta(id_n)); 
						FEg( id ) += gK[m][n](z,x) * real(Dv_eta_dag(id_m) * Du_eta(id_n)); 
						FEg( id ) += gK[m][n](x,z) * real(Du_eta_dag(id_m) * Dv_eta(id_n)); 
					}
				}
				FE_minus_ref += wx*wz*(h*h * (FEb(id) - FEdens(id)) + 2.0/3.0*FEg( id )); 
			}
		}

		FEg *= 2.0/3.0/(h*h);
		FEdens = FEb + FEg; 
		if(Nu==1 || Nv==1) FE_minus_ref /=h;

		return FE_minus_ref;
	}
// ===========================================================


// ===========================================================
// FiveCompHe3 :: function definitions
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

	double FiveCompHe3::FreeEn(T_vector & OPvector, in_conditions cond, Eigen::VectorXd & FEdens, Eigen::VectorXd & FEb, Eigen::VectorXd & FEg) {
		T_vector dummy(vect_size); 
		this->bulkRHS_FE(cond, OPvector, dummy, FEb); // get the bulk contribution to free energy density

		T_vector Du_eta(vect_size), Dv_eta(vect_size); 
		if( Nu > 1 ) Du_eta = Du_FE * OPvector; else Du_eta = T_vector::Zero(vect_size);
		if( Nv > 1 ) Dv_eta = Dv_FE * OPvector; else Dv_eta = T_vector::Zero(vect_size);

		T_vector Du_eta_dag=Du_eta.adjoint(), Dv_eta_dag=Dv_eta.adjoint(); 
		
		double FE_minus_ref=0.0; // relative to the reference Free energy (FEdens on INPUT) 
		int x = 0, z = 1;            // indices for the K-matrix
		double wx, wz;               // weights for the integral free energy (trapezoidal rule)

		for (int u = 0; u < Nu; u++) {
			if( (u==0 || u==Nu-1) && Nu>1 ) wx=0.5; else wx=1.0;
			for (int v = 0; v < Nv; v++) {
				if( (v==0 || v==Nv-1) && Nv>1 ) wz=0.5; else wz=1.0;

				int id = ID(u,v,0);
				FEg( id ) = 0;

				for (int m = 0; m < Nop; m++) {
					for (int n = 0; n < Nop; n++) {
						int id_m = ID(u,v,m);
						int id_n = ID(u,v,n);

						FEg( id ) += gK[m][n](x,x) * real(Du_eta_dag(id_m) * Du_eta(id_n)); 
						FEg( id ) += gK[m][n](z,z) * real(Dv_eta_dag(id_m) * Dv_eta(id_n)); 
						FEg( id ) += gK[m][n](z,x) * real(Dv_eta_dag(id_m) * Du_eta(id_n)); 
						FEg( id ) += gK[m][n](x,z) * real(Du_eta_dag(id_m) * Dv_eta(id_n)); 
					}
				}
				FE_minus_ref += wx*wz*(h*h * (FEb(id) -FEdens(id)) + 2.0/3.0*FEg( id )); 
			}
		}

		FEg *= 2.0/3.0/(h*h);
		FEdens = FEb + FEg; 
		if(Nu==1 || Nv==1) FE_minus_ref /=h;

		return FE_minus_ref;
	}

	void FiveCompHe3 :: initGuessWithCircularDomain(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, vector<int> & no_update) {
		// use "conditions5_AyyFlip.txt"

		// all boundaries should be Dirichlet => 1
		// we will create the circular domain in the middle of Ayy to be -1

		// the radius of the domain...we're assuming that the domain itself is large enough
		double r = 20;
		if (Nu*h/2 <= r || Nv*h/2 <= r)
			cout << "WARNING: the circular domain will be too large!" << endl;
		
		// location of the ceter of the grid
		int u_center = round(Nu/2);
		int v_center = round(Nv/2);

		// loop over the whole mesh
		for (int u = 0; u < Nu; u++) { 
			double x= h*(u-u_center); 
			for (int v = 0; v < Nv; v++) { 
				double z= h*(v-v_center);
				double rho = sqrt(x*x + z*z);

				OPvector(ID(u,v,0))=1;  
				OPvector(ID(u,v,1))=tanh( (rho -r)/2.0 ); 
				OPvector(ID(u,v,2))=1;
				OPvector(ID(u,v,3))=0;
				OPvector(ID(u,v,4))=0;

				// no update for the Dirichlet boundaries and OP1 at the radius of the circle domain
				for(int n=0; n<5; n++){ 
					int id = ID(u,v,n);
					if(u==0  && eta_BC[n].typeL=="D" && Nu>1) no_update.push_back( id );
					if(u==Nu-1 && eta_BC[n].typeR=="D" && Nu>1) no_update.push_back( id );
					if(v==0  && eta_BC[n].typeB=="D" && Nv>1) no_update.push_back( id );
					if(v==Nv-1 && eta_BC[n].typeT=="D" && Nv>1) no_update.push_back( id );
					if( n==1 &&  abs(rho-r) < 0.4*h) no_update.push_back( id );
				}
			}
		}
		return;
	}
// ===========================================================


// ===========================================================
// OneCompSC :: function definitions
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

	double OneCompSC::FreeEn(T_vector & OPvector, in_conditions cond, Eigen::VectorXd & FEdens, Eigen::VectorXd & FEb, Eigen::VectorXd & FEg) {
		T_vector dummy(vect_size); 
		this->bulkRHS_FE(cond, OPvector, dummy, FEb); // get the bulk contribution to free energy density

		T_vector Du_eta(vect_size), Dv_eta(vect_size); 
		if( Nu > 1 ) Du_eta = Du_FE * OPvector; else Du_eta = T_vector::Zero(vect_size);
		if( Nv > 1 ) Dv_eta = Dv_FE * OPvector; else Dv_eta = T_vector::Zero(vect_size);

		T_vector Du_eta_dag=Du_eta.adjoint(), Dv_eta_dag=Dv_eta.adjoint(); 
		
		double FE_minus_ref=0.0; // relative to a reference OP configuration 
		int x = 0, z = 1;            // indices for the K-matrix
		double wx, wz;               // weights for the integral free energy (trapezoidal rule)

		for (int u = 0; u < Nu; u++) {
			if( (u==0 || u==Nu-1) && Nu>1 ) wx=0.5; else wx=1.0;
			for (int v = 0; v < Nv; v++) {
				if( (v==0 || v==Nv-1) && Nv>1 ) wz=0.5; else wz=1.0;

				int id = ID(u,v,0);
				FEg( id ) = 0;

				FEg( id ) += gK[0][0](x,x) * real(Du_eta_dag(id) * Du_eta(id)); 
				FEg( id ) += gK[0][0](z,z) * real(Dv_eta_dag(id) * Dv_eta(id)); 
				FEg( id ) += gK[0][0](z,x) * real(Dv_eta_dag(id) * Du_eta(id)); 
				FEg( id ) += gK[0][0](x,z) * real(Du_eta_dag(id) * Dv_eta(id)); 

				FE_minus_ref += wx*wz*(h*h * (FEb(id) - FEdens(id)) + 2.0*FEg( id )); 
			}
		}

		FEg *= 2.0/(h*h);
		FEdens = FEb + FEg; 
		if(Nu==1 || Nv==1) FE_minus_ref /=h;

		return FE_minus_ref;
	}
// ===========================================================


// ===========================================================
// Cylindrical :: function definitions
	void Cylindrical::Build_curvilinear_matrices(std::vector<Bound_Cond> eta_BC) {
		Dr.resize(9*grid_size,9*grid_size);
		Dz.resize(9*grid_size,9*grid_size);
		Dphi.resize(9*grid_size,9*grid_size);
		r_inv.resize(grid_size,grid_size); // make sure to only pre-multiply by r_inv!

		// calcualte the 1/r - matrix		
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			int id = ID(u,v,0);
			r_inv.coeffRef(id,id) = 1./(u+u_shift); 
		}

		T_vector dummy; // just a dummy variable for arguments below
		SpMat_cd tempDv, tempDu;
		
		for (int n = 0; n < 9; n++) {
			if (Nv > 1) {
				tempDv = n < Nop ? Dv_BD(eta_BC[n],n,dummy,dummy) : Dv; // we only have BC's for the Nop ...
				Dz += Place_subMatrix(n,n,9,tempDv);
			}
			if (Nu > 1) {
				tempDu = n < Nop ? Du_BD(eta_BC[n],n,dummy,dummy) : Du; // ... so we just use the default otherwise
				Dr += Place_subMatrix(n,n,9,tempDu);
			}
		}
		
		// build the D_phi matrix
		MatrixXd temp(9,9);
		temp << 	0, 0, 0, 0, 0,-1,-1, 0, 0,
				0, 0, 0, 0, 0, 1, 1, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 0, 0,-1,
				0, 0, 0, 0, 0, 0, 0,-1, 0,
				1,-1, 0, 0, 0, 0, 0, 0, 0,
				1,-1, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 1, 0, 0, 0, 0,
				0, 0, 0, 1, 0, 0, 0, 0, 0;
		Dphi = kroneckerProduct(temp,r_inv);
		
		return;
	}

	void Cylindrical::BuildSolverMatrixCyl( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, std::vector<Bound_Cond> eta_BC) {
		// make all the matrices in eq. (54)
		SpMat_cd Dr2t(grid_size,grid_size), // D_r^2 ~
				Drp(grid_size,grid_size),   // D_r^+
				Drm(grid_size,grid_size),   // D_r^-
				Drzp(grid_size,grid_size),  // D_rz^+
				D2t(grid_size,grid_size);   // D^2 ~
				// initialize to 0-matrix in case they aren't needed
				//    to be built, but they will still be in the code!
				Dr2t *= 0.0;
				Drp *= 0.0;
				Drm *= 0.0;
				Drzp *= 0.0;
				D2t *= 0.0;

		// define D matrices that are used in others' definitions
		SpMat_cd Dr2_gsize(grid_size,grid_size),
				Dr_gsize(grid_size,grid_size),
				Dr_over_r(grid_size,grid_size),
				Drz_gsize(grid_size,grid_size),
				Dz_over_r(grid_size,grid_size),
				Dz2_gsize(grid_size,grid_size),
				Dz_gsize(grid_size,grid_size);
				Dr2_gsize *= 0.0;
				Dr_gsize *= 0.0;
				Dr_over_r *= 0.0;
				Drz_gsize *= 0.0;
				Dz_over_r *= 0.0;
				Dz2_gsize *= 0.0;
				Dz_gsize *= 0.0;

		SpMat_cd r2_inv = r_inv*r_inv; // this one is the same regardless of BC's

		// the matrices -- 2 parts of the solver matrix
		SpMat_cd DK23(vect_size,vect_size), DK1(vect_size,vect_size);

		// a temporary rhsBC to add to the actual
		T_vector rhsBC_r2_gsize = T_vector::Zero(vect_size);
		T_vector rhsBC_z2_gsize = T_vector::Zero(vect_size);
		T_vector rhsBC_rz_gsize = T_vector::Zero(vect_size);
		T_vector rhsBC_r_gsize = T_vector::Zero(vect_size);
		T_vector rhsBC_z_gsize = T_vector::Zero(vect_size);

		for (int n = 0; n < Nop; n++) {
			// initialize the D matrices based on BC's
			// only make those that we ABOLUTELY need, because these can take a LONG time to build!
			if (Nu > 1) Dr2_gsize = Du2_BD(eta_BC[n],n,initOPvector,rhsBC_r2_gsize);
			if (Nv > 1) Dz2_gsize = Dv2_BD(eta_BC[n],n,initOPvector,rhsBC_z2_gsize);
			if (Nu > 1 && Nv > 1 && (n == 0 || n == 2 || n == 3 || n == 4)) Drz_gsize = Duv_BD(eta_BC[n],n,initOPvector,rhsBC_rz_gsize);
			if (Nu > 1 && (n == 0 || n == 1)) Dr_gsize  = Du_BD (eta_BC[n],n,initOPvector,rhsBC_r_gsize);
			if (Nv > 1 && n == 0) Dz_gsize  = Dv_BD (eta_BC[n],n,initOPvector,rhsBC_z_gsize);
			if (Nu > 1) Dr_over_r = r_inv*Dr_gsize;
			if (Nv > 1 && (n == 0 || n == 1 || n == 3 || n == 4)) Dz_over_r = r_inv*Dz_gsize;

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
				DK23 += Place_subMatrix(0,n,Nop,Dr2t)
					+ Place_subMatrix(1,n,Nop,r_inv*Drp);
				if (Nop > 3) DK23 += Place_subMatrix(4,n,Nop,Drzp);
				
				DK1 += Place_subMatrix(0,n,Nop,D2t-r2_inv)
					+ Place_subMatrix(1,n,Nop,2.*r2_inv);
				
			} else if (n == 1) {
				DK23 += Place_subMatrix(0,n,Nop,-r_inv*Drm)
					+ Place_subMatrix(1,n,Nop,-r2_inv);
				if (Nop > 3) DK23 += Place_subMatrix(4,n,Nop,-Dz_over_r);
				
				DK1 += Place_subMatrix(1,n,Nop,D2t-r2_inv)
					+ Place_subMatrix(0,n,Nop,2.*r2_inv);

			} else if (n == 2) {
				DK23 += Place_subMatrix(2,n,Nop,Dz2_gsize);
				if (Nop > 3) DK23 += Place_subMatrix(3,n,Nop,Drz_gsize);
				
				DK1 += Place_subMatrix(n,n,Nop,D2t+r2_inv);

			} else if (n == 3) {
				DK23 += Place_subMatrix(2,n,Nop,Drzp)
					+ Place_subMatrix(3,n,Nop,Dr2t);
				
				DK1 += Place_subMatrix(n,n,Nop,D2t);

			} else if (n == 4) {
				DK23 += Place_subMatrix(0,n,Nop,Drz_gsize)
					+ Place_subMatrix(1,n,Nop,Dz_over_r)
					+ Place_subMatrix(4,n,Nop,Dz2_gsize);
				
				DK1 += Place_subMatrix(n,n,Nop,D2t);
			}
		}

		// make M from the 2 matrices
		M = K1*DK1 + K23*DK23; // use the K-values defined in "Kstorage.hpp"

		return;
	}

	// this one should be the same as the 5-component system, right? (says just below eq. (38).)
	void Cylindrical::bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & FEb) {
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

	double Cylindrical::FreeEn(T_vector & OPvector, in_conditions parameters, Eigen::VectorXd & FEdensity, Eigen::VectorXd & freeEb, Eigen::VectorXd & freeEg) {
		T_vector dummy(vect_size);
		this->bulkRHS_FE(parameters, OPvector, dummy, freeEb);

		// We MUST have eta be the full size so that we can most easily account for
		//	 the four D_phi contributions from components we don't calculate. (If we
		//   don't make it the full size, then we will either get mis-matched matrix
		//   multiplication or we will miss the "extra" parts that come from D_phi.)
		T_vector eta(9*grid_size); eta << OPvector, T_vector::Zero((9-Nop)*grid_size);
		T_vector eta_c=eta.conjugate();
		
		double FE_minus_ref=0.0; // relative to a reference Free energy, FEdensity on INPUT 
		double wr, wz;               // weights for the integral free energy (trapezoidal rule)
		// some other variable used in the loops below; label them here to reduce time in the loops
		vector<T_vector> D_eta_c, D_eta;
		double K_ijmn, r;
		int id, id_m, id_n;

		// build the derivative vectors 
		D_eta_c.insert(D_eta_c.end(),{Dr*eta_c, Dphi*eta_c, Dz*eta_c});
		D_eta.insert(D_eta.end(),{Dr*eta, Dphi*eta, Dz*eta});

		for (int u = 0; u < Nu; u++) { // loop over the whole grid
			if( Nu>1 && ( u==Nu-1 || (u==0 && u_shift<0.1) ) ) wr=0.5; else wr=1.0;
			r = (u+u_shift); 
			for (int v = 0; v < Nv; v++) {
				if( (v==0 || v==Nv-1) && Nv>1 ) wz=0.5; else wz=1.0;

				id = ID(u,v,0);
				freeEg(id)=0;
				for (int m = 0; m < 9; m++) for (int n = 0; n < 9; n++) { // sum over op components 
					id_m = ID(u,v,m), id_n = ID(u,v,n); 
					for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { // sums over coordinates 
						K_ijmn = Kij(i,j,m,n); 
						if ( abs(K_ijmn) > 1e-4 )   freeEg(id) += K_ijmn * D_eta_c[i](id_m) * D_eta[j](id_n); // see eq. (46)
				}}

				FE_minus_ref += (2.0*M_PI*r*h) * wr*wz * ( h*h*(freeEb(id) - FEdensity(id)) + 2.0/3.0*freeEg(id) ); // see eq. (47)

			} // v
		} // u

		freeEg *= 2.0/3.0/(h*h);
		FEdensity = freeEb + freeEg;
		if(Nu==1 || Nv==1) FE_minus_ref /=h; 

		return FE_minus_ref;
	}

	// initial guess functions
	void Cylindrical::initialOPguess_Cylindrical_bubble(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, std::vector<int> & no_update) {
		// build the guess from the 3 component cylindrical solution
		int Nop_init = 3;
		string file_name_init = "conditions3c.txt";
		in_conditions cond_init;
		vector<Bound_Cond> eta_BC_init;
		read_input_data(Nop_init, cond_init, eta_BC_init, file_name_init);
		cond_init.SIZEu = Nu;
		cond_init.SIZEv = Nv;
		cond_init.STEP = h;
		vector<int> no_update_init;
		int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
		int VectSize_init = cond_init.Nop * GridSize_init;
		T_vector OPvector_init(VectSize_init);
		T_vector rhsBC_init = T_vector::Zero(VectSize_init);
		SpMat_cd M_init(VectSize_init,VectSize_init);
		SC_class *pSC_init;
		pSC_init = new Cylindrical( Nop_init, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP, eta_BC_init );
		pSC_init->initialOPguess_Cylindrical_simple3(eta_BC_init, OPvector_init, no_update_init);
		pSC_init->BuildSolverMatrixCyl( M_init, rhsBC_init, OPvector_init, eta_BC_init );
		Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init); // get the solution to the 3-comp. system

		OPvector *= 0.0; // zero out the OPvector

		// go through the first 2 OP components, copying values from the solution
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			int id1 = ID(u,v,0);
			int id2 = ID(u,v,1);
			OPvector(id1) = OPvector_init(id1);
			OPvector(id2) = OPvector_init(id2);
		}

		double radius, r_wall = 0.6*min(Nu,Nv)*h; // may need to change r_wall...; r0, z0, 
		// going through the entire grid
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			double r = h*(u+u_shift); // current r coordinate value
			double z = h*v;           // current z value
			int id = ID(u,v,2);       // just look at the Azz component right now
			radius = sqrt(r*r + z*z); // radius of the bubble
			OPvector( id ) = tanh( (radius-r_wall)/5. ) * tanh( z/5. );
		}

		// // set some boundary values based on conditions file
		// for (int n = 0; n < Nop; n++)
		// if  (n != 2)
		// for (int u = 0; u < Nu; u++) // loop over the whole mesh
		// for (int v = 0; v < Nv; v++) {
		// 	int id = ID(u,v,n);

		// 	if (u == 0 && eta_BC[n].typeL == string("D")) {
		// 		OPvector(id) = eta_BC[n].valueL;
		// 		no_update.push_back(id);
		// 	} else if (u == Nu-1 && eta_BC[n].typeR == string("D")) {
		// 		OPvector(id) = eta_BC[n].valueR;
		// 		no_update.push_back(id);
		// 	} else if (v == 0 && eta_BC[n].typeB == string("D")) {
		// 		OPvector(id) = eta_BC[n].valueB;
		// 		no_update.push_back(id);
		// 	} else if (v == Nv-1 && eta_BC[n].typeT == string("D")) {
		// 		OPvector(id) = eta_BC[n].valueT;
		// 		no_update.push_back(id);	
		// 	}
		// 	else {//if (n < 2)
		// 		// build a smooth guess based on BC's
		// 		int wT = v, wL = Nu-u, wB = Nv-v, wR = u; // weights
		// 		OPvector(id) = ( eta_BC[n].valueB * wB
		// 						+ eta_BC[n].valueT * wT
		// 						+ eta_BC[n].valueL * wL
		// 						+ eta_BC[n].valueR * wR )
		// 					/(Nu+Nv);
		// 	}
		// 	// else OPvector(id) = 0.0;
		// }

		VectorXd dummyFE=VectorXd::Zero(grid_size); 
		this->WriteAllToFile(OPvector, dummyFE, dummyFE, dummyFE, "initial_guess_OP5c.txt");
		delete pSC_init;

		return;
	}

	// just the simplest system: only z-variation, no domain walls, just the pairbreaking at the surface
	void Cylindrical::initialOPguess_Cylindrical_simple3(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, std::vector<int> & no_update) {
		// Using BC's to make a smooth guess
		for (int n = 0; n < Nop; n++)
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			int id = ID(u,v,n);
			if (n < 2) {
				// build a smooth guess based on BC's
				int wT = v, wL = Nu-u, wB = Nv-v, wR = u; // weights
				OPvector(id) = ( eta_BC[n].valueB * wB
								+ eta_BC[n].valueT * wT
								+ eta_BC[n].valueL * wL
								+ eta_BC[n].valueR * wR )
							/(Nu+Nv);
			} else { // n == 2
				double z = h*v; // so that z = 0 is the surface
				OPvector(id) = tanh(z/5.0);
			}

			if (u == 0 && eta_BC[n].typeL == string("D")) {
				OPvector(id) = eta_BC[n].valueL;
				no_update.push_back(id);
			}
			else if (u == Nu-1 && eta_BC[n].typeR == string("D")) {
				OPvector(id) = eta_BC[n].valueR;
				no_update.push_back(id);
			}
			else if (v == 0 && eta_BC[n].typeB == string("D")) {
				OPvector(id) = eta_BC[n].valueB;
				no_update.push_back(id);
			}
			else if (v == Nv-1 && eta_BC[n].typeT == string("D")) {
				OPvector(id) = eta_BC[n].valueT;
				no_update.push_back(id);	
			}
		}

		VectorXd dummyFE(grid_size); 
		this->WriteAllToFile(OPvector, dummyFE, dummyFE, dummyFE, "initial_guess_OP3c.txt");
	}

	void Cylindrical::initialOPguess_Cylindrical_simple5(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, std::vector<int> & no_update) {
		for (int n = 0; n < Nop; n++)
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			int id = ID(u,v,n);
			if (n < 2) {
				// build a smooth guess based on BC's
				int wT = v, wL = Nu-u, wB = Nv-v, wR = u; // weights
				OPvector(id) = ( eta_BC[n].valueB * wB
								+ eta_BC[n].valueT * wT
								+ eta_BC[n].valueL * wL
								+ eta_BC[n].valueR * wR )
							/(Nu+Nv);
			} else if (n == 2) {
				double z = h*v; // so that z = 0 is the surface
				OPvector(id) = tanh(z/5.0);
			} else { // n == 3, 4
				OPvector(id) = 0.0;
			}

			if (u == 0 && eta_BC[n].typeL == string("D")) {
				OPvector(id) = eta_BC[n].valueL;
				no_update.push_back(id);
			}
			else if (u == Nu-1 && eta_BC[n].typeR == string("D")) {
				OPvector(id) = eta_BC[n].valueR;
				no_update.push_back(id);
			}
			else if (v == 0 && eta_BC[n].typeB == string("D")) {
				OPvector(id) = eta_BC[n].valueB;
				no_update.push_back(id);
			}
			else if (v == Nv-1 && eta_BC[n].typeT == string("D")) {
				OPvector(id) = eta_BC[n].valueT;
				no_update.push_back(id);	
			}
		}

		VectorXd dummyFE(grid_size); 
		this->WriteAllToFile(OPvector, dummyFE, dummyFE, dummyFE, "initial_guess_OP5c.txt");
	}

	void Cylindrical::initialOPguess_Cylindrical_AzzFlip(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, std::vector<int> & no_update) {
		for (int n = 0; n < Nop; n++)
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			int id = ID(u,v,n);
			if (n < 2) {
				// build a smooth guess based on BC's
				// int wT = v, wL = Nu-u-1, wB = Nv-v-1, wR = u; // weights
				OPvector(id) = 1.0;
							// 	( eta_BC[n].valueB * wB
							// 	+ eta_BC[n].valueT * wT
							// 	+ eta_BC[n].valueL * wL
							// 	+ eta_BC[n].valueR * wR )
							// /(Nu+Nv);
			} else if (n == 2) {
				double r = h*(u+u_shift);
				double r_center = h*Nu/2;
				OPvector(id) = tanh((r-r_center)/5.0);
			} else { // n == 3, 4
				OPvector(id) = 0.0;
			}

			if (u == 0 && eta_BC[n].typeL == string("D")) {
				OPvector(id) = eta_BC[n].valueL;
				no_update.push_back(id);
			}
			else if (u == Nu-1 && eta_BC[n].typeR == string("D")) {
				OPvector(id) = eta_BC[n].valueR;
				no_update.push_back(id);
			}
			else if (v == 0 && eta_BC[n].typeB == string("D")) {
				OPvector(id) = eta_BC[n].valueB;
				no_update.push_back(id);
			}
			else if (v == Nv-1 && eta_BC[n].typeT == string("D")) {
				OPvector(id) = eta_BC[n].valueT;
				no_update.push_back(id);	
			}
		}

		VectorXd dummyFE(grid_size); 
		this->WriteAllToFile(OPvector, dummyFE, dummyFE, dummyFE, "initial_guess_OP5c.txt");
	}
// ===========================================================
