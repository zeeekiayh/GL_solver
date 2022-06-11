#include "SC_classes.hpp"
#include "readwrite.hpp"
#include "structures.hpp"
#include <vector>
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

	SpMat_cd Place_subMatrix(int i, int j, int size, SpMat_cd sm);
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

	void ThreeCompHe3::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[]) 
	{
		T_vector Du_eta(vect_size), Dv_eta(vect_size); 

		if( Nu > 1 ) Du_eta = Du_FE * OPvector; else Du_eta = T_vector::Zero(vect_size);
		if( Nv > 1 ) Dv_eta = Dv_FE * OPvector; else Dv_eta = T_vector::Zero(vect_size);

		T_vector Du_eta_dag=Du_eta.adjoint(), 
				Dv_eta_dag=Dv_eta.adjoint(); 
		
		int x = 0, z = 1; // indeces for the K-matrix

		for (int u = 0; u < Nu; u++) {
			for (int v = 0; v < Nv; v++) {
				int FE_id = ID(u,v,0);
				freeEg( FE_id ) = 0;

				for (int m = 0; m < Nop; m++) {
					for (int n = 0; n < Nop; n++) {
						int id_m = ID(u,v,m);
						int id_n = ID(u,v,n);

						freeEg( FE_id ) += gK[m][n](x,x) * (Du_eta_dag(id_m) * Du_eta(id_n)).real(); 
						freeEg( FE_id ) += gK[m][n](z,z) * (Dv_eta_dag(id_m) * Dv_eta(id_n)).real(); 
						freeEg( FE_id ) += gK[m][n](z,x) * (Dv_eta_dag(id_m) * Du_eta(id_n)).real(); 
						freeEg( FE_id ) += gK[m][n](x,z) * (Du_eta_dag(id_m) * Dv_eta(id_n)).real(); 
					}
				}
			}
		}
		freeEg *= 2.0/3.0/(h*h);
		return; 
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

	void FiveCompHe3::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[]) 
	{
		T_vector Du_eta(vect_size), Dv_eta(vect_size); 

		if( Nu > 1 ) Du_eta = Du_FE * OPvector; else Du_eta = T_vector::Zero(vect_size);
		if( Nv > 1 ) Dv_eta = Dv_FE * OPvector; else Dv_eta = T_vector::Zero(vect_size);

		T_vector Du_eta_dag=Du_eta.adjoint(), 
				Dv_eta_dag=Dv_eta.adjoint(); 
		
		int x = 0, z = 1; // indeces for the K-matrix

		for (int u = 0; u < Nu; u++) {
			for (int v = 0; v < Nv; v++) {
				int FE_id = ID(u,v,0);
				freeEg( FE_id ) = 0;

				for (int m = 0; m < Nop; m++) {
					for (int n = 0; n < Nop; n++) {
						int id_m = ID(u,v,m);
						int id_n = ID(u,v,n);

						freeEg( FE_id ) += gK[m][n](x,x) * (Du_eta_dag(id_m) * Du_eta(id_n)).real(); 
						freeEg( FE_id ) += gK[m][n](z,z) * (Dv_eta_dag(id_m) * Dv_eta(id_n)).real(); 
						freeEg( FE_id ) += gK[m][n](z,x) * (Dv_eta_dag(id_m) * Du_eta(id_n)).real(); 
						freeEg( FE_id ) += gK[m][n](x,z) * (Du_eta_dag(id_m) * Dv_eta(id_n)).real(); 
					}
				}
			}
		}
		freeEg *= 2.0/3.0/(h*h);
		return; 
	}

	double FiveCompHe3 :: defectEnergy(const Eigen::VectorXd & freeEb, const Eigen::VectorXd & freeEg) {
		double DefectEnergy = 0.;
		for (int u = 0; u < Nu; u++) {
			for (int v = 0; v < Nv; v++) {
				int FE_id = ID(u,v,0);
				DefectEnergy += freeEg(FE_id);
			}
		}
		return DefectEnergy * (h*h);
	}

	void FiveCompHe3 :: initGuessWithCircularDomain(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
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

	void OneCompSC::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[]) {

		T_vector Du_eta(vect_size), Dv_eta(vect_size); 

		if( Nu > 1 ) Du_eta = Du_FE * OPvector; else Du_eta = T_vector::Zero(vect_size);
		if( Nv > 1 ) Dv_eta = Dv_FE * OPvector; else Dv_eta = T_vector::Zero(vect_size);

		T_vector Du_eta_dag=Du_eta.adjoint(), 
				Dv_eta_dag=Dv_eta.adjoint(); 
		
		int x = 0, z = 1; // indeces for the K-matrix

		for (int u = 0; u < Nu; u++) {
			for (int v = 0; v < Nv; v++) {
				int FE_id = ID(u,v,0);
				freeEg( FE_id ) = 0;

					int id_m = ID(u,v,0);
					int id_n = ID(u,v,0);

					freeEg( FE_id ) += gK[0][0](x,x) * (Du_eta_dag(id_m) * Du_eta(id_n)).real(); 
					freeEg( FE_id ) += gK[0][0](z,z) * (Dv_eta_dag(id_m) * Dv_eta(id_n)).real(); 
					freeEg( FE_id ) += gK[0][0](z,x) * (Dv_eta_dag(id_m) * Du_eta(id_n)).real(); 
					freeEg( FE_id ) += gK[0][0](x,z) * (Du_eta_dag(id_m) * Dv_eta(id_n)).real(); 
			}
		}
		freeEg *= 2/(h*h);
		return; 
	}
// ===========================================================


//* ===========================================================
// Cylindrical :: function definitions
	void Cylindrical::Build_curvilinear_matrices(Bound_Cond eta_BC[]) {
		// cout << "\tin Build_curvilinear_matrices" << endl;
		Ident.resize(vect_size,vect_size); // set the identity matrix to the right size
		Ident.setIdentity();
		
		// cout << "\there 0" << endl;

		Dr.resize(vect_size,vect_size);
		Dz.resize(vect_size,vect_size);
		Dphi.resize(vect_size,vect_size);
		r_inv.resize(grid_size,grid_size);

		// cout << "\there 1" << endl;

		T_vector initOPvector; // ? should it actually be passed in as an argument in this function?
		T_vector rhsBC; // ? should it actually be passed in as an argument?
		// neither of these are actually used below...we just have them as arguments in
		// 	the D*_BD() functions below, but those functions don't use them either!
		
		// cout << "\there 2" << endl;

		// SpMat_cd R(grid_size,grid_size); // there is no need for the vect_size x vect_size r_inv?
		for (int u = 0; u < Nu; u++) {
			for (int v = 0; v < Nv; v++) {
				int id = ID(u,v,0);
				r_inv.coeffRef(id,id) = (u == 0) ? 1e10 : 1./(u*h);
			}
		}
		
		// cout << "\there 3" << endl;

		for (int n = 0; n < Nop; n++) {
			// cout << "\t\tn=" << n << endl;
			Dr += Place_subMatrix(n,n,Nop,Du_BD(eta_BC[n],n,initOPvector,rhsBC));
			Dz += Place_subMatrix(n,n,Nop,Dv_BD(eta_BC[n],n,initOPvector,rhsBC));
			// r_inv += Place_subMatrix(n,n,Nop,R);
		}
		
		// cout << "\there 4" << endl;

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
		
		// cout << "\there 5" << endl;

		return;
	}

	void Cylindrical::BuildSolverMatrixCyl( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[]) {
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
				DK23 += Place_subMatrix(0,n,Nop,Dr2t)
					+ Place_subMatrix(1,n,Nop,Drp*r_inv)
					+ Place_subMatrix(4,n,Nop,Drzp);
				
				DK1 += Place_subMatrix(0,n,Nop,D2t-r2_inv)
					+ Place_subMatrix(1,n,Nop,2.*r2_inv);
			} else if (n == 1) {
				DK23 += Place_subMatrix(0,n,Nop,-Drm)
					+ Place_subMatrix(1,n,Nop,-r2_inv)
					+ Place_subMatrix(4,n,Nop,-Dz_over_r);
				
				DK1 += Place_subMatrix(1,n,Nop,D2t-r2_inv)
					+ Place_subMatrix(0,n,Nop,2.*r2_inv);
			} else if (n == 2) {
				DK23 += Place_subMatrix(2,n,Nop,Dz2_gsize)
					+ Place_subMatrix(3,n,Nop,Drz_gsize);
				
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

		// make M from these 2 matrices
		M = 1.*DK1 + 2.*DK23; // TODO: FIX K-VALUES HERE!
		// Should they be passed in as arguments?

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

	void Cylindrical::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[]) {
		// we only need this for the energy inspection...not for the solution
	}

	double Cylindrical::defectEnergy(const Eigen::VectorXd & freeEb, const Eigen::VectorXd & freeEg) {
		//
		return 0.;
	}

	// initial guess functions
	void Cylindrical::initialOPguess_Cylindrical_bubble(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update) {
		// Nop, Nu, Nv, h, ID() - are part of the SC_class, so we use them!
		double radius, r_wall = 0.5*min(Nu,Nv)*h; // may need to change r_wall...; r0, z0, 

		// for (int n = 0; n < Nop; n++) {
			// going through the entire grid
			for (int u = 0; u < Nu; u++) {
				double r = h*u; 

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

	void Cylindrical::initialOPguess_Cylindrical_simple(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update) {
		// TODO
	}

	void Cylindrical::initialOPguess_Cylindrical_AzzFlip(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update) {
		// TODO
	}
// ===========================================================*/