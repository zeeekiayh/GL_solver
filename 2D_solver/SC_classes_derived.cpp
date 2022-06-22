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

	double ThreeCompHe3::FreeEn(T_vector & OPvector, in_conditions cond, Eigen::VectorXd & FEdens, Eigen::VectorXd & FEb, Eigen::VectorXd & FEg) {
		T_vector dummy(vect_size); 
		this->bulkRHS_FE(cond, OPvector, dummy, FEb); // get the bulk contribution to free energy density

		T_vector Du_eta(vect_size), Dv_eta(vect_size); 
		if( Nu > 1 ) Du_eta = Du_FE * OPvector; else Du_eta = T_vector::Zero(vect_size);
		if( Nv > 1 ) Dv_eta = Dv_FE * OPvector; else Dv_eta = T_vector::Zero(vect_size);

		T_vector Du_eta_dag=Du_eta.adjoint(), Dv_eta_dag=Dv_eta.adjoint(); 
		
		double FE_minus_uniform=0.0; // relative to uniform state with density FEuniform=-1;
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
				FE_minus_uniform += wx*wz*(h*h * (FEb(id) + 1.0) + 2.0/3.0*FEg( id )); 
			}
		}

		FEg *= 2.0/3.0/(h*h);
		FEdens = FEb + FEg; 
		if(Nu==1 || Nv==1) FE_minus_uniform /=h;

		return FE_minus_uniform;
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
		
		double FE_minus_uniform=0.0; // relative to uniform state with density FEuniform=-1;
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
				FE_minus_uniform += wx*wz*(h*h * (FEb(id) + 1.0) + 2.0/3.0*FEg( id )); 
			}
		}

		FEg *= 2.0/3.0/(h*h);
		FEdens = FEb + FEg; 
		if(Nu==1 || Nv==1) FE_minus_uniform /=h;

		return FE_minus_uniform;
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
		
		double FE_minus_uniform=0.0; // relative to uniform state with density FEuniform=-1;
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

				FE_minus_uniform += wx*wz*(h*h * (FEb(id) + 1.0) + 2.0*FEg( id )); 
			}
		}

		FEg *= 2.0/(h*h);
		FEdens = FEb + FEg; 
		if(Nu==1 || Nv==1) FE_minus_uniform /=h;

		return FE_minus_uniform;
	}
// ===========================================================


// ===========================================================
// Cylindrical :: function definitions
	void Cylindrical::Build_curvilinear_matrices(std::vector<Bound_Cond> eta_BC) {
		Dr.resize(vect_size,vect_size);
		Dz.resize(vect_size,vect_size);
		Dphi.resize(vect_size,vect_size);

		// make sure to only pre-multiply by r_inv!
		r_inv.resize(grid_size,grid_size);
		r_inv_full.resize(vect_size,vect_size);

		T_vector initOPvector; // ? should it actually be passed in as an argument in this function?
		T_vector rhsBC; // ? should it actually be passed in as an argument?
		// neither of these are actually used below...we just have them as arguments in
		// 	the D*_BD() functions below, but those functions don't use them either!
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			int id = ID(u,v,0);
			r_inv.coeffRef(id,id) = 1./((u+0.5)*h); // shift u-points to half-grid points to avoid 1/r|r=0 problems
		}
		
		for (int n = 0; n < Nop; n++) {
			            Dz += Place_subMatrix(n,n,Nop,Dv_BD(eta_BC[n],n,initOPvector,rhsBC));
			if (Nu > 1) Dr += Place_subMatrix(n,n,Nop,Du_BD(eta_BC[n],n,initOPvector,rhsBC));
			r_inv_full += Place_subMatrix(n,n,Nop,r_inv); // we do need this one for the free energy calculations...see eq.'s (44 - 47)
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
		
		return;
	}

	void Cylindrical::BuildSolverMatrixCyl( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, std::vector<Bound_Cond> eta_BC) {
		cout << "BuildSolverMatrixCyl" << endl;
		// make all the matrices in eq. (54)
		SpMat_cd Dr2t(grid_size,grid_size), // D_r^2 ~
				Drp(grid_size,grid_size),   // D_r^+
				Drm(grid_size,grid_size),   // D_r^-
				Drzp(grid_size,grid_size),  // D_rz^+
				D2t(grid_size,grid_size);   // D^2 ~

		// define D matrices that are used in others' definitions
		SpMat_cd Dr2_gsize(grid_size,grid_size),
				Dr_gsize(grid_size,grid_size),
				Dr_over_r(grid_size,grid_size),
				Drz_gsize(grid_size,grid_size),
				Dz_over_r(grid_size,grid_size),
				Dz2_gsize(grid_size,grid_size),
				Dz_gsize(grid_size,grid_size);
		SpMat_cd r2_inv = r_inv*r_inv; // this one is the same regardless of BC's

		// test
		// Dr_gsize *= 0.0;

		// the matrices -- 2 parts of the solver matrix
		SpMat_cd DK23(vect_size,vect_size), DK1(vect_size,vect_size);

		// a temporary rhsBC to add to the actual
		T_vector rhsBC_r2_gsize = T_vector::Zero(vect_size);
		T_vector rhsBC_z2_gsize = T_vector::Zero(vect_size);
		T_vector rhsBC_rz_gsize = T_vector::Zero(vect_size);
		T_vector rhsBC_r_gsize = T_vector::Zero(vect_size);
		T_vector rhsBC_z_gsize = T_vector::Zero(vect_size);

		// Should they be passed in as arguments?
		double K1_const = 1.0, K23_const = 2.0; // WHAT TO DO WITH K-VALUES?!

		for (int n = 0; n < Nop; n++) {
			// initialize the D matrices based on BC's
			// only make those that we ABOLUTELY need, because these take a LONG time to build!
			if (Nu > 1) Dr2_gsize = Du2_BD(eta_BC[n],n,initOPvector,rhsBC_r2_gsize);
			// Dr2_gsize *= 0.0;
			Dz2_gsize = Dv2_BD(eta_BC[n],n,initOPvector,rhsBC_z2_gsize);
			if (Nu > 1 && (n == 0 || n == 2 || n == 3 || n == 4)) Drz_gsize = Duv_BD(eta_BC[n],n,initOPvector,rhsBC_rz_gsize);
			if (Nu > 1 && (n == 0 || n == 1)) Dr_gsize  = Du_BD (eta_BC[n],n,initOPvector,rhsBC_r_gsize);
			if (n == 0) Dz_gsize  = Dv_BD (eta_BC[n],n,initOPvector,rhsBC_z_gsize);
			if (Nu > 1) Dr_over_r = r_inv*Dr_gsize;
			if (n == 0 || n == 1 || n == 3 || n == 4) Dz_over_r = r_inv*Dz_gsize;

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
				// rhsBC += K23_const*rhsBC_r2_gsize;
				// rhsBC += K23_const*rhsBC_r_gsize;

				if (Nop > 3) {
					DK23 += Place_subMatrix(4,n,Nop,Drzp);
					// rhsBC += K23_const*rhsBC_rz_gsize;
				}
				
				DK1 += Place_subMatrix(0,n,Nop,D2t-r2_inv)
					+ Place_subMatrix(1,n,Nop,2.*r2_inv);
				// rhsBC += K1_const*rhsBC_z2_gsize;
				// rhsBC += K1_const*rhsBC_r2_gsize;
			} else if (n == 1) {
				DK23 += Place_subMatrix(0,n,Nop,-r_inv*Drm)
					+ Place_subMatrix(1,n,Nop,-r2_inv);
				// rhsBC += K23_const*rhsBC_r_gsize;

				if (Nop > 3)
					DK23 += Place_subMatrix(4,n,Nop,-Dz_over_r);
				
				DK1 += Place_subMatrix(1,n,Nop,D2t-r2_inv)
					+ Place_subMatrix(0,n,Nop,2.*r2_inv);
				// rhsBC += K1_const*rhsBC_z2_gsize;
				// rhsBC += K1_const*rhsBC_r2_gsize;
			} else if (n == 2) {
				DK23 += Place_subMatrix(2,n,Nop,Dz2_gsize);
				// rhsBC += K23_const*rhsBC_z2_gsize;
				if (Nop > 3) {
					DK23 += Place_subMatrix(3,n,Nop,Drz_gsize);
					// rhsBC += K23_const*rhsBC_rz_gsize;
				}
				
				DK1 += Place_subMatrix(n,n,Nop,D2t+r2_inv);
				// rhsBC += K1_const*rhsBC_z2_gsize;
				// rhsBC += K1_const*rhsBC_r2_gsize;
			} else if (n == 3) {
				DK23 += Place_subMatrix(2,n,Nop,Drzp)
					+ Place_subMatrix(3,n,Nop,Dr2t);
				// rhsBC += K23_const*rhsBC_rz_gsize;
				// rhsBC += K23_const*rhsBC_r2_gsize;
				
				DK1 += Place_subMatrix(n,n,Nop,D2t);
				// rhsBC += K1_const*rhsBC_z2_gsize;
				// rhsBC += K1_const*rhsBC_r2_gsize;
			} else if (n == 4) {
				DK23 += Place_subMatrix(0,n,Nop,Drz_gsize)
					+ Place_subMatrix(1,n,Nop,Dz_over_r)
					+ Place_subMatrix(4,n,Nop,Dz2_gsize);
				// rhsBC += K23_const*rhsBC_rz_gsize;
				// rhsBC += K23_const*rhsBC_z2_gsize;
				
				DK1 += Place_subMatrix(n,n,Nop,D2t);
				// rhsBC += K1_const*rhsBC_z2_gsize;
				// rhsBC += K1_const*rhsBC_r2_gsize;
			}

		}

		// make M from the 2 matrices
		M = K1_const*DK1 + K23_const*DK23;

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
		// TODO
		return 0.0;
	}

	// initial guess functions
	void Cylindrical::initialOPguess_Cylindrical_bubble(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, std::vector<int> & no_update) {
		// Nop, Nu, Nv, h, ID() - are part of the SC_class, so we use them!
		double radius, r_wall = 0.5*min(Nu,Nv)*h; // may need to change r_wall...; r0, z0, 

		// for (int n = 0; n < Nop; n++) {
			// going through the entire grid
			for (int u = 0; u < Nu; u++) {
				double r = h*(u+0.5);  // shift u-points to half-grid points // old: double r = h*u;

				for (int v = 0; v < Nv; v++) {
					double z = h*v;  

					int id = ID(u,v,2);
					radius = sqrt(r*r + z*z);
					OPvector( id ) = tanh( (radius-r_wall)/2. );
				}
			}
		// }

		// set some boundary values based on conditions file
		for (int n = 0; n < Nop; n++)
		if  (n != 2)
		for (int u = 0; u < Nu; u++) // loop over the whole mesh
		for (int v = 0; v < Nv; v++) {
			int id = ID(u,v,n);

			if (u == 0) {
				OPvector(id) = eta_BC[n].valueL;
			} else if (u == Nu-1) {
				OPvector(id) = eta_BC[n].valueR;
			} else if (v == 0) {
				OPvector(id) = eta_BC[n].valueB;
			} else if (v == Nv-1) {
				OPvector(id) = eta_BC[n].valueT;
			}
		}
		// smooth off the initial guess a little
		for (int i = 0; i < 10; i++)
		for (int n = 0; n < Nop; n++)
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			if ( n != 2 ) {
				int id = ID(u,v,n);
				int idU = ID(u,(v<Nv-1)?v+1:v,n);
				int idL = ID((u>0)?u-1:u,v,n);
				int idD = ID(u,(v>0)?v-1:v,n);
				int idR = ID((u<Nu-1)?u+1:u,v,n);
				OPvector(id) = (OPvector(idU) + OPvector(idL) + OPvector(idD) + OPvector(idR))/4.;
			}
		}

		return;
	}

	// just the simplest system: only z-variation, no domain walls, just the pairbreaking at the surface
	void Cylindrical::initialOPguess_Cylindrical_simple3(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, std::vector<int> & no_update) {
		cout << "initialOPguess_Cylindrical_simple3" << endl;
		// start with the simple 3-comp. system in cartesian
		// SC_class *pSC;
		// int Nop3 = 3;
		// pSC = new ThreeCompHe3(Nop3,Nu,Nv,h);

		// in_conditions cond3;
		// Bound_Cond *BC3;
		// read_input_data(Nop3,cond3,BC3,"conditions3_normal.txt");
		// vector<int> no_update3;
		// int GridSize = Nu * Nv;
		// int VectSize = Nop3 * GridSize;
		// T_vector OPvector3(VectSize);
		// T_vector rhsBC = T_vector::Zero(VectSize);
		// SpMat_cd M(VectSize,VectSize);

		// cond3.SIZEu = Nu;
		// cond3.SIZEv = Nv;
		// cond3.STEP = h;

		// pSC->initialOPguess(BC3,OPvector3,no_update3);
		// pSC->BuildSolverMatrix( M, rhsBC, OPvector3, BC3 );
		// Solver(OPvector3, M, rhsBC, cond3, no_update3, pSC); // solve the system setup above
		// OPvector = OPvector3;

		// Using BC's to make a smooth guess
		for (int n = 0; n < Nop; n++)
		for (int u = 0; u < Nu; u++)
		for (int v = 0; v < Nv; v++) {
			int id = ID(u,v,n);
			if (n < 2) {
				// OPvector(id) = 1.0;
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
			} //else {
				//OPvector(id) = OPvector3(id);
			//}
		}

		VectorXd freeEb(grid_size), freeEg(grid_size);
		this->WriteAllToFile(OPvector, freeEb, freeEg, "initial_guess_OP3.txt");
	}

	void Cylindrical::initialOPguess_Cylindrical_simple5(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, std::vector<int> & no_update) {
		// start with the simple 3-comp. system in cartesian
		SC_class *pSC;
		int Nop3 = 3;
		pSC = new ThreeCompHe3(Nop3,Nu,Nv,h);

		in_conditions cond3;
		vector<Bound_Cond> BC3;
		read_input_data(Nop3,cond3,BC3,"conditions3_normal.txt");
		vector<int> no_update3;
		int GridSize = Nu * Nv;
		int VectSize = Nop3 * GridSize;
		T_vector OPvector3(VectSize);
		T_vector rhsBC = T_vector::Zero(VectSize);
		SpMat_cd M(VectSize,VectSize);

		cond3.SIZEu = Nu;
		cond3.SIZEv = Nv;
		cond3.STEP = h;

		pSC->initialOPguess(BC3,OPvector3,no_update3);
		pSC->BuildSolverMatrix( M, rhsBC, OPvector3, BC3 );
		Solver(OPvector3, M, rhsBC, cond3, no_update3, pSC); // solve the system setup above
		
		OPvector = T_vector::Zero(vect_size); // size for the 5 comp
		for (int id = 0; id < VectSize; id++) OPvector(id) = OPvector3(id);

		// for (int n = 0; n < Nop; n++)
		// for (int u = 0; u < Nu; u++)
		// for (int v = 0; v < Nv; v++) {
		// 	int id = ID(u,v,n);
		// 	if (n < 2) {
		// 		// build a smooth guess based on BC's
		// 		int wT = v, wL = Nu-u, wB = Nv-v, wR = u; // weights
		// 		OPvector(id) = ( eta_BC[n].valueB * wB
		// 						+ eta_BC[n].valueT * wT
		// 						+ eta_BC[n].valueL * wL
		// 						+ eta_BC[n].valueR * wR )
		// 					/(Nu+Nv);
		// 	} else if (n == 2) {
		// 		double z = h*v; // so that z = 0 is the surface
		// 		OPvector(id) = tanh(z/5.0);
		// 	} else
		// 		OPvector(id) = 0.0;

		// 	if (u == 0 && eta_BC[n].typeL == string("D"))
		// 		no_update.push_back(id);
		// 	else if (u == Nu-1 && eta_BC[n].typeR == string("D"))
		// 		no_update.push_back(id);
		// 	else if (v == 0 && eta_BC[n].typeB == string("D"))
		// 		no_update.push_back(id);
		// 	else if (v == Nv-1 && eta_BC[n].typeT == string("D"))
		// 		no_update.push_back(id);
		// }
	}

	void Cylindrical::initialOPguess_Cylindrical_AzzFlip(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, std::vector<int> & no_update) {
		// TODO
	}
// ===========================================================