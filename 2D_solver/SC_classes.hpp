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
		Eigen::Matrix2d **gK;            // gradient coefficients in GL functional, each gK matrix is (xx, xz, // zx, zz)
      	SpMat_cd Du2, Dv2, Duv, Du, Dv; // base derivative matrices, without boundary conditions
		SpMat_cd Du_FE, Dv_FE; // full first-order Derivative matrices to be used in FE calculation 
	public:
		//constructor 
		// SC_op () {Nop=1; *eta_bulk = 1;}; // default constructor for single-component OP 
		//      - eliminate default constructor for now, and always initiate the class ourselves  
		SC_class(int n, int nx, int ny, double step) {
			Nop=n; Nu=nx; Nv=ny; grid_size=nx*ny; vect_size=grid_size * Nop; h=step; 
			eta_bulk=new T_scalar [n];
			gK = new Eigen::Matrix2d *[n]; 
			for (int i = 0; i < n; i++) gK[i] = new Eigen::Matrix2d [n];
			this->Build_D_Matrices(); // build all general D-matrices for anything that will need them
		};

		//destructor: declared  virtual so that derived classes could destroy the base class
		virtual ~SC_class () {
			delete[] eta_bulk;
			for(int i = 0; i <Nop; i++) delete[] gK[i]; 
			delete[] gK; 
		}

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
		void initOPguess_special(in_conditions cond, Bound_Cond eta_BC[], T_vector & OPvector, Eigen::Matrix2d **gradK, std::vector<int> & no_update); 
		//void initGuessWithCircularDomain(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update);
		// not sure why we need virtual here but OK ... -Anton 
		// Needs to be virtual so that it can be accessed by "pSC" in 'gl_fdm.cpp', because it is defined as a 'SC_class' object;
		// 		any initial guess functions that should be defined for specific OP's will need to them be defined in the derived
		//		classes, and thus will all need to be virtual as well.
		virtual void initGuessWithCircularDomain(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update){};
		virtual void initialOPguess_Cylindrical (Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update){};

		// the general method of building the solver matrix
		void BuildSolverMatrix   ( SpMat_cd & M, T_vector & rhsBC, const T_vector initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK );
		virtual void BuildSolverMatrixCyl( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK){};
		
		// write out a vector to a file
		void WriteToFile(const T_vector& vector, std::string file_name, int flag);
		void WriteAllToFile(const T_vector& solution, const T_vector& FE_bulk, const T_vector& FE_grad, std::string file_name);

		// to be defined by each derived class
		virtual void bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & FEb){};
		virtual void gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK){};
		virtual void gradFE_curv(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK){};
		virtual double defectEnergy(const Eigen::VectorXd & freeEb, const Eigen::VectorXd & freeEg){return 0.;};
};

// derived_classes depends on the particular form of the bulk free energy, that gives the appropriate RHS vector
class OneCompSC : public SC_class {
	protected:
		// some new variables... 
	public: 
		//constructor: we need to construct the base class too 
		OneCompSC (int n, int nx, int ny, double step) : SC_class(n, nx, ny, step) {
				gK[0][0] << 1, 0, 
						0, 1; 
		} 
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
		ThreeCompHe3 (int n, int nx, int ny, double step) : SC_class(n, nx, ny, step) {
			int c[]={0,1,2}; // fill the gK matrix for OP components 0,1,2=(Axx, Ayy, Azz)
			for(int m=0; m<3; m++) for(int k=0; k<3; k++){
				gK[m][k] << Kuu[c[m]][c[k]], Kuv[c[m]][c[k]], 
						Kvu[c[m]][c[k]], Kvv[c[m]][c[k]] ;
				//std::cout << "\n gK["<<m<<"]["<<k<<"]=\n" << gK[m][k] << "\n";
			}
		} 
		~ThreeCompHe3 () {} 

		// void BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK); 
		void bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & freeEb);
		void gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK);

		// iniital guess function prototypes
		void initOPguess_1DNarrowChannel(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update);
};

class FiveCompHe3 : public SC_class {
	protected:
	public: 
		//constructor: we need to construct the base class too 
		FiveCompHe3 (int n, int nx, int ny, double step) : SC_class(n, nx, ny, step) { 
			int c[]={0,1,2,3,4}; // fill the gK matrix for OP components 0,1,2,3,4=(Axx, Ayy, Azz, Azx, Axz)
			for(int m=0; m<5; m++) for(int k=0; k<5; k++){
				gK[m][k] << Kuu[c[m]][c[k]], Kuv[c[m]][c[k]], 
						Kvu[c[m]][c[k]], Kvv[c[m]][c[k]] ;
			}
		}
		~FiveCompHe3 () {} 

		// void BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK); 
		void bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & freeEb);
		void gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK);
		double defectEnergy(const Eigen::VectorXd & freeEb, const Eigen::VectorXd & freeEg);

		// initial guess function prototypes
		void initGuessWithCircularDomain(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update);
};

class Cylindrical : public SC_class {
	protected:
		// inherits all the basic D matrices
		
		// define covariant matrices ::: I can be used for imaginary unit sometimes 
		SpMat_cd Dr, Dz, Dphi, Ident, r_inv;
	public:
		// the constructor should call the parent constructor and some other functions.
		Cylindrical (int n, int nr, int nz, double h, Bound_Cond eta_BC[]) : SC_class(n,nr,nz,h) {
			// std::cout << "In Cylynd constructor" << std::endl;
			// call derivative-matrix-building functions here
			Build_curvilinear_matrices(eta_BC);
			// call K-matrix-building functions here? --No, we do that elsewhere
		}
		~Cylindrical(){}

		void Build_curvilinear_matrices(Bound_Cond eta_BC[]);
		// void BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK); 
		void BuildSolverMatrixCyl( SpMat_cd & M, T_vector & rhsBC, const T_vector & initOPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK);
		void bulkRHS_FE(in_conditions parameters, T_vector & OPvector, T_vector & newRHSvector, Eigen::VectorXd & freeEb);
		void gradFE(Eigen::VectorXd & freeEg, const T_vector & OPvector, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK);
		double defectEnergy(const Eigen::VectorXd & freeEb, const Eigen::VectorXd & freeEg);

		// initial guess function prototypes
		void initialOPguess_Cylindrical(Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update);
};

#endif
