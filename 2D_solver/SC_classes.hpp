#ifndef SC_class_hpp_
#define SC_class_hpp_

// the typedefs and custom structures are in this file:
#include "structures.hpp"

// base_class of superconducting OP; contains the common features 
class SC_class{
	protected: // protected variables can be accessed by the derived classes
		int Nop; // number of OP components
		int Nu, Nv, grid_size; // size of the grid in u,v-directions, and total number of grid points 
		int vect_size; // full vector size of 
		double h; // grid step size
		T_scalar *eta_bulk; // OP components (pointer) 
      	SpMat_cd Du2, Dv2, Duv, Du, Dv; // base derivative matrices, without boundary conditions
      	// bool use_no_update = 0; // says whether or not to use the no_update vector
      	// vector<int> no_update; // stores all the indices that will not be modified in the RHS
		SpMat_cd FEgrad; // the F_hat_grad matrix in equation 36 in Latex file
	public:
		//constructor 
		// SC_op () {Nop=1; *eta_bulk = 1;}; // default constructor for single-component OP 
		//      - eliminate default constructor for now, and always initiate the class ourselves  
		SC_class(int n, int nx, int ny, double step) {
			Nop=n; Nu=nx; Nv=ny; grid_size=nx*ny; vect_size=grid_size * Nop; h=step; 
			eta_bulk=new T_scalar [n];
			this->Build_D_Matrices(); // build all general D-matrices for anything that will need them
		}; 

		//destructor: declared  virtual so that derived classes could destroy the base class
		virtual ~SC_class () {delete eta_bulk;}

		// common functions -----------------------------------------------------------
		int sizeOP (void) {return Nop;}

		// unique id for each grid point (n_u,n_v) and order parameter component n
		// we assign ID by taking OP component and then going lift-right, and down-up; repeat for next OP component  
		int ID(int n_u, int n_v, int n) { return n_u + Nu*n_v + grid_size*n; } 

		// building the general derivative matrices except boundary points 
		void Build_D_Matrices();

		// functions that build derivative matrices that incorporate boundary conditions 
		// -- will be called from derived class functions 
		SpMat_cd Du2_BD(Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC);
		SpMat_cd Dv2_BD(Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC);
		SpMat_cd Duv_BD(Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC);
		SpMat_cd  Du_BD(Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC);
		SpMat_cd  Dv_BD(Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC);

		// the general method of making the initial guess based on the given BC's
		void initialOPguess(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update);
		// a method of initializing a guess based on a previous solution
		// void initialOPguessFromSolution(SC_class & SC, Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update);
		void initialOPguessFromSolution(const T_vector & solution, T_vector & OPvector, std::vector<int> & no_update);

		// the general method of building the solver matrix
		void BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK );
		
		// write out a vector to a file
		void WriteToFile(const T_vector& vector, std::string file_name, int flag);

		// to be defined by each derived class
		virtual void bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & FEb){};
		virtual void gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK){};
};

// derived_classes depends on the particular form of the bulk free energy, that gives the appropriate RHS vector
class OneCompSC : public SC_class {
	protected:
		// some new variables... 
	public: 
		//constructor: we need to construct the base class too 
		OneCompSC (int n, int nx, int ny, double step) : SC_class(n, nx, ny, step) {} 
		//destructor : automatically calls virtual ~SC_class()
		~OneCompSC () {} 

		// void BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK); 
		void bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & freeEb);
		void gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK);
};

class ThreeCompHe3 : public SC_class {
	protected:
	public: 
		//constructor: we need to construct the base class too 
		ThreeCompHe3 (int n, int nx, int ny, double step) : SC_class(n, nx, ny, step) {} 
		~ThreeCompHe3 () {} 

		// void BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK); 
		void bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & freeEb);
		void gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK);
};

class FiveCompHe3 : public SC_class {
	protected:
	public: 
		//constructor: we need to construct the base class too 
		FiveCompHe3 (int n, int nx, int ny, double step) : SC_class(n, nx, ny, step) {} 
		~FiveCompHe3 () {} 

		// void BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK); 
		void bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & freeEb);
		void gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK);
};

#endif
