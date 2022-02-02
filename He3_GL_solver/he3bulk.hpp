#ifndef _he3bulk
#define _he3bulk

#include <cmath>
#include <complex>
#include <eigen/Eigen/Dense>

using std::complex;
using namespace Eigen;

// bulk values of He3 free energy and of dFdA 3x3 matrix 
// normalized by the bulk values of the He3-B gap and the free energy 
// at temperature t and pressure p
void he3bulk(double, double, Matrix3cd, Matrix3cd&, double&);

// beta_i/beta_B parameters 
// beta_B  = beta12 + beta345/3
// include (t,p)-dependent strong coupling effects here 
void betas(double, double, double*);

// base_class of superconducting OP; contains the common features 
class SC_op{
	protected:
		int Nop; // number of OP components
	public:
		complex<double> *eta_bulk; // pointer to bulk OP structure 

		//constructor 
		SC_op (); // default constructor for single-component OP
		SC_op (int, complex<double>*); 
		//destructor
		/*virtual*/ ~SC_op () {delete eta_bulk;}

		// common functions 
		int size (void);
		// some functions for gradient terms that take K-matrix and Bound conditions as input
		// ..... 

		// functions that get replaced by appropriate functions of derived classes ("virtual")
		virtual double bulkRHS(double, double, complex<double>*, complex<double>*, double&); 
};

// derived_class depends on the particular form of the bulk free energy, that gives the appropriate RHS vector
class FiveCompHe3 : public SC_op {
	protected:
		// some new variables 
	public: 
		//constructor: we need to construct the base class too 
		FiveCompHe3 (int, complex<double>*);
		// --- need to add default constructor ... SC_op () {}; ????
		//destructor : automatically calls ~SC_op() ?
		~FiveCompHe3();

		double bulkRHS(double, double, complex<double>*, complex<double>*, double&);
};

#endif