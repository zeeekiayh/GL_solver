#include "SC_classes.hpp"
#include "readwrite.hpp"
#include "structures.hpp"
#include <vector>
#include <eigen/unsupported/Eigen/KroneckerProduct> // for the Kronecker product in Place_subMatrix

using namespace std;
using namespace Eigen;

// ===========================================================
// Function prototypes defined in other files
	// prototype for function defined in linear_eq_solver.cpp
	void Solver(T_vector & f, SpMat_cd M, T_vector rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC);

	// function to get the RHS of GL differential equations, determined from the bulk functional; 
	// ------  coded in he3bulk.cpp ---------
	// input: temp t, pressure p, matrix A, 
	// output: betaB parameter that gives the bulk free energy of B-phase; derivative matrix dFdA, bulk free energy FE
	void he3bulk(double t, double p, double & betaB, Eigen::Matrix<T_scalar,3,3> A, Eigen::Matrix<T_scalar,3,3> & dFdA, double & FE);

	// prototype for function defined in 'SC_classes_base.cpp'
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

	void ThreeCompHe3::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
		// Follow equation 10 and 36 in Latex doc
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
	}

	void ThreeCompHe3 :: initOPguess_1DNarrowChannel(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
		// use "conditions3_1DChannel.txt"
		// cout << "initOPguess_NarrowChannel()" << endl;
		// solve for the normal 3-component system
		int Nop_init = 3;
		in_conditions cond_init;
		Bound_Cond eta_BC_init[Nop_init];
		Matrix2d **gradK_init;
		gradK_init = new Matrix2d *[Nop_init];
		for (int i = 0; i < Nop_init; i++) gradK_init[i] = new Matrix2d [Nop_init];
		read_input_data(Nop_init, cond_init, eta_BC_init, "conditions3_normal.txt");

		// gradK_init

		cond_init.maxStore = 5;
		cond_init.rel_p = 0.1;
		cond_init.wait = 1;
		vector<int> no_update_init;
		// for us, we need this to be the z-length of our system
		cond_init.SIZEv = Nv;
		cond_init.SIZEu = 1;
		cond_init.STEP = h;
		int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
		int VectSize_init = cond_init.Nop * GridSize_init;
		T_vector OPvector_init(VectSize_init);
		SpMat_cd M_init(VectSize_init,VectSize_init);
		T_vector rhsBC_init = T_vector::Zero(VectSize_init);
		SC_class *pSC_init = new ThreeCompHe3( Nop_init, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP );
		pSC_init->initialOPguess(eta_BC_init, OPvector_init, no_update_init);
		pSC_init->BuildSolverMatrix( M_init, rhsBC_init, OPvector_init, eta_BC_init, gradK_init );
		Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init);

		for (int n = 0; n < Nop; n++) {
			for (int v = 0; v < Nv; v++) {
				int id_forward = ID(0,v,n);
				int id_backward = ID(0,Nv-1-v,n);
				OPvector(id_forward) = OPvector_init(id_forward)*OPvector_init(id_backward);
			}
		}
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

	// a method of initializing a guess based on a previous solution
	// NOTE: this funciton will mess things up if grid_size for solution is different than the grid_size for OPvector
	void FiveCompHe3 :: initialOPguessFromSolution(T_vector & OPvector, std::vector<int> & no_update) {
		cout << "initialOPguessFromSolution" << endl;
		// To calculate a solution based on the initial guess of this here...
		// -----------------------------------------------------------------------------
		int Nop_init = 3;
		in_conditions cond_init;
		Bound_Cond eta_BC_init[Nop_init];
		// Matrix2d **gradK_init;
		// gradK_init = new Matrix2d *[Nop_init];
		// for (int i = 0; i < Nop_init; i++) gradK_init[i] = new Matrix2d [Nop_init];

		read_input_data(Nop_init, cond_init, eta_BC_init, "conditions"+to_string(Nop_init)+".txt");

		vector<bool> components;
		components.insert(components.end(),{true,false,false,true}); // this may need to be changed for your system!
		auto gradK_init = kMatrix::smallKMatrix(Nop_init, components);

		// confirm_input_data(Nop_init, cond_init, eta_BC_init, gradK_init);

		cond_init.maxStore = 5;
		cond_init.rel_p = 0.1;
		cond_init.wait = 1;

		vector<int> no_update_init;

		// make sure they are the right sizes!
		cond_init.SIZEv = Nv;
		cond_init.SIZEu = Nu;
		cond_init.STEP = h;
		int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
		int VectSize_init = cond_init.Nop * GridSize_init;

		T_vector OPvector_init(VectSize_init);
		T_vector rhsBC_init = T_vector::Zero(VectSize_init);
		SpMat_cd M_init(VectSize_init,VectSize_init);

		SC_class *pSC_init;
		pSC_init = new ThreeCompHe3( Nop_init, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP );

		pSC_init->initialOPguess(eta_BC_init, OPvector_init, no_update_init);
		pSC_init->BuildSolverMatrix( M_init, rhsBC_init, OPvector_init, eta_BC_init, gradK_init );
		Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init);
			// -----------------------------------------------------------------------------

		if (OPvector_init.size() > OPvector.size()) {
			cout << "ERROR: can't initialize a guess with a previous solution larger than the current:\n\tsolution.size() = " << OPvector_init.size() << "; OPvector.size()" << OPvector.size() << endl;
			return;
		}
		// should we fix the boundaries...? using no_update?
		for (int i = 0; i < OPvector_init.size(); i++)
			OPvector(i) = OPvector_init(i);
		for (int i = OPvector_init.size(); i < OPvector.size(); i++)
			OPvector(i) = 0.;
		// And what will happen to the solution of the new one (i.e. what OPvector will become)
		//    if the derivative matrices are different? Will it still converge?
	}

	void FiveCompHe3 :: initOPguess_AzzFlip_WS2016(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
		// use "conditions5_W&S2016.txt"
		// solve for the normal 3-component system
		int Nop_init = 3;
		in_conditions cond_init;
		Bound_Cond eta_BC_init[Nop_init];
		Matrix2d **gradK_init;
		gradK_init = new Matrix2d *[Nop_init];
		for (int i = 0; i < Nop_init; i++) gradK_init[i] = new Matrix2d [Nop_init];
		read_input_data(Nop_init, cond_init, eta_BC_init, "conditions3_1DChannel.txt");

		// gradK_init

		cond_init.maxStore = 5;
		cond_init.rel_p = 0.1;
		cond_init.wait = 1;
		vector<int> no_update_init;
		// for us, we need this to be the z-length of our system
		cond_init.SIZEv = Nv;
		cond_init.SIZEu = 1;
		cond_init.STEP = h;
		int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
		int VectSize_init = cond_init.Nop * GridSize_init;
		T_vector OPvector_init(VectSize_init);
		SpMat_cd M_init(VectSize_init,VectSize_init);
		T_vector rhsBC_init = T_vector::Zero(VectSize_init);
		SC_class *pSC_init = new ThreeCompHe3( Nop_init, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP );
			pSC_init->initOPguess_1DNarrowChannel(eta_BC_init, OPvector_init, no_update_init);
		pSC_init->BuildSolverMatrix( M_init, rhsBC_init, OPvector_init, eta_BC_init, gradK_init );
		Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init);

		// going through the entire grid
		for (int u = 0; u < Nu; u++) {
			for (int v = 0; v < Nv; v++) {
				for (int n = 0; n < Nop; n++) {
					double x = h*u;
						int id      = ID(u,v,n);
					int id_init = v + GridSize_init*n;

					if (n == 0 || n == 1)
					OPvector(id) = OPvector_init( id_init );
					else if (n == 2)
					OPvector(id) = tanh((x-h*Nu/2)/2) * OPvector_init( id_init );
					else if (n == 3 || n == 4)
					OPvector(id) = 0.;
				}
			}
		}
		return;
	}

	// "and the interesting thing to check is take the initial guess of 3-compnent solution on the
	//    left side smoothly changing to 3-component solution with Azz flipped on the right side."
	void FiveCompHe3 :: initOPguess_AzzFlip(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
		// use "conditions5_xDeform.txt"
		// solve for the normal 3-component system for the initial guess
		int Nop_init = 3;
		in_conditions cond_init;
		Bound_Cond eta_BC_init[Nop_init];
		Matrix2d **gradK_init;
		gradK_init = new Matrix2d *[Nop_init];
		for (int i = 0; i < Nop_init; i++) gradK_init[i] = new Matrix2d [Nop_init];
		read_input_data(Nop_init, cond_init, eta_BC_init, "conditions3_normal.txt");

		// gradK_init
		
		cond_init.maxStore = 5;
		cond_init.rel_p = 0.1;
		cond_init.wait = 1;
		vector<int> no_update_init;
		// for us, we need this to be the z-length of our system
		cond_init.SIZEv = Nv;
		cond_init.SIZEu = 1;
		cond_init.STEP = h;
		int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
		int VectSize_init = cond_init.Nop * GridSize_init;
		T_vector OPvector_init(VectSize_init);
		SpMat_cd M_init(VectSize_init,VectSize_init);
		T_vector rhsBC_init = T_vector::Zero(VectSize_init);
		SC_class *pSC_init = new ThreeCompHe3( Nop_init, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP );
		pSC_init->initialOPguess(eta_BC_init, OPvector_init, no_update_init);
		pSC_init->BuildSolverMatrix( M_init, rhsBC_init, OPvector_init, eta_BC_init, gradK_init );
		Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init);

		for (int n = 0; n < Nop; n++) {
			for (int v = 0; v < Nv; v++) {
				double z = h*v/16.;
				for (int u = 0; u < Nu; u++) {
					double x = h*(u-Nu/2)/4.;
					int id = ID(u,v,n);
					int id_init = v + GridSize_init*n;
					
					if (n < 3) {
						OPvector(id) = OPvector_init(id_init) * ( (n==2)?tanh(x):1. );
					}// else 
					if (n==2||n == 3 || n == 4) {
						if (n == 4) {
							OPvector(id) = 0.5 * pow(2.718, -x/2.*x/2.)*-x/2. * z*pow(2.718, -z*z);
						} else if (n == 3) {
							OPvector(id) = 0.;
						}

						if (u == 0) {
							OPvector(id) = eta_BC[n].valueL;
							if (eta_BC[n].typeL == "D") { no_update.push_back(id); /*cout << "adding id=" << id << " to 'no_update'." << endl;*/ }
						} else if (u == Nu-1) {
							OPvector(id) = eta_BC[n].valueR;
							if (eta_BC[n].typeR == "D") { no_update.push_back(id); /*cout << "adding id=" << id << " to 'no_update'." << endl;*/ }
						} else if (v == 0) {
							OPvector(id) = eta_BC[n].valueB;
							if (eta_BC[n].typeB == "D") { no_update.push_back(id); /*cout << "adding id=" << id << " to 'no_update'." << endl;*/ }
						} else if (v == Nv-1) {
							OPvector(id) = eta_BC[n].valueT;
							if (eta_BC[n].typeT == "D") { no_update.push_back(id); /*cout << "adding id=" << id << " to 'no_update'." << endl;*/ }
						}
					}
				}
			}
		}
		return;
	}

	void FiveCompHe3 :: initGuessWithCircularDomain(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
		// use "conditions5_AyyFlip.txt"

		// all boundaries should be Dirichlet => 1
		// we will create the circular domain in the middle of Ayy to be -1

		// the radius of the domain...we're assuming that the domain itself is large enough
		double r = 16;
		if (Nu*h/2 <= r || Nv*h/2 <= r)
			cout << "WARNING: the circular domain will be too large!" << endl;
		
		// location of the ceter of the grid
		int u_center = round(Nu/2);
		int v_center = round(Nv/2);

		for (int n = 0; n < Nop; n++) {
			// loop over the whole mesh
			for (int u = 0; u < Nu; u++) {
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
					} else if (n == 1) { // for only Ayy:
						OPvector(id) = tanh(  (sqrt(pow(h*(u-u_center),2) + pow(h*(v-v_center),2)) - r)/3.  );
					} else if (n == 3 || n == 4) {
						OPvector(id) = 0.;
					} else {
						OPvector(id) = 1.;
					}
				}
			}
		}
		// // smooth off the initial guess a little
		// for (int i = 0; i < 30; i++) {
		// 	for (int n = 0; n < Nop; n++) {
		// 	   for (int u = 0; u < Nu; u++) {
		// 	      for (int v = 0; v < Nv; v++) {
		// 	         if ( (u > 0) && (u < Nu-1) && (v > 0) && (v < Nv-1) && (n == 0) && (n == 2) ) {
		// 	            int id = ID(u,v,n);
		// 	               // cout << "id = " << id << endl;
		// 	            int idU = ID(u,v+1,n);
		// 	               // cout << "idU = " << idU << endl;
		// 	            int idL = ID(u-1,v,n);
		// 	               // cout << "idL = " << idL << endl;
		// 	            int idD = ID(u,v-1,n);
		// 	               // cout << "idD = " << idD << endl;
		// 	            int idR = ID(u+1,v,n);
		// 	               // cout << "idR = " << idR << endl;
		// 	            OPvector(id) = (OPvector(idU) + OPvector(idL) + OPvector(idD) + OPvector(idR))/4.;
		// 	            // cout << "OPvector(id) = " << (OPvector(idU) + OPvector(idL) + OPvector(idD) + OPvector(idR))/4. << endl;
		// 	         }
		// 	      }
		// 	   }
		// 	}
		// }
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
	void OneCompSC::gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) {
		//
	}
// ===========================================================


// ===========================================================
// Cylindrical :: function definitions
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
		return 0.;
	}
// ===========================================================