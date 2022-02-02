#include "he3bulk.hpp"

// Function to get the RHS of GL differential equations, determined from the bulk functional;
// Bulk values of He3 free energy and of dFdA 3x3 matrix 
//  normalized by the bulk values of the He3-B gap and the free energy 
//  at temperature t and pressure p
void he3bulk(double t, double p, Matrix3cd A, Matrix3cd & dFdA, double & FE)
{
	double beta[6];
	auto Adag = A.adjoint(); // Hermitian conjugate
	auto At = A.transpose(); // transpose 
	auto Ac = A.conjugate(); // complex conjugate

	betas(t,p, beta);

	// derivative of FE w.r.to A^* -components 
	dFdA = beta[0]*A; 
	dFdA 	+= 1.0/3.0*beta[1] * (A*At).trace() *Ac;
	dFdA 	+= 1.0/3.0*beta[2] * (A*Adag).trace() *A;
	dFdA 	+= 1.0/3.0*beta[3] * (A * At * Ac);
	dFdA 	+= 1.0/3.0*beta[4] * (A * Adag * A);
	dFdA 	+= 1.0/3.0*beta[5] * (Ac * At * A);

	// bulk free energy, normalized by B-phase value
	FE = 2.0/3.0 * beta[0] * real( (A * Adag).trace() );
	FE 	+= 1.0/9.0 * beta[1] * norm(  (A*At).trace()  );
	FE 	+= 1.0/9.0 * beta[2] * norm( (A*Adag).trace() );
	FE 	+= 1.0/9.0 * beta[3] * real( (A * At * Ac * Adag ).trace() );
	FE 	+= 1.0/9.0 * beta[4] * real( (A * Adag * A * Adag).trace() );
	FE 	+= 1.0/9.0 * beta[5] * real( (A * Adag * Ac * At ).trace() );

	return;
}

// beta_i/beta_B parameters 
// beta_B  = beta12 + beta345/3
// include (t,p)-dependent strong coupling effects here 
void betas(double t, double p, double *beta)
{
	beta[0]=-1.0; // alpha-coefficient in GL / |alpha|
	// strong-coupling corrections will change beta[1-5]
	beta[1]=1.0; // reference beta, followed by weak-coupling relations
	beta[2]=-2*beta[1];
	beta[3]=-2*beta[1];
	beta[4]=-2*beta[1];
	beta[5]=+2*beta[1];

	double betaB = beta[1]+beta[2] + (beta[3]+beta[4]+beta[5])/3.0;

	for(int i=1; i<=5; i++) beta[i] /= betaB;
	
	return;
}



// base_class of superconducting OP; contains the common features
//constructor 
SC_op::SC_op () {Nop=1; *eta_bulk = 1;}; // default constructor for single-component OP
SC_op::SC_op (int OP_size, complex<double> *etab) {}
//destructor
// virtual SC_op::~SC_op () {delete eta_bulk;}

// common functions 
int SC_op::size (void) {return Nop;}
// some functions for gradient terms that take K-matrix and Bound conditions as input
// ..... 

// functions that get replaced by appropriate functions of derived classes ("virtual")
// virtual void SC_op::bulkRHS(double t, double p, complex<double> *eta, complex<double> *dFdeta, double & FEbulk);



// derived_class depends on the particular form of the bulk free energy, that gives the appropriate RHS vector
//constructor: we need to construct the base class too 
FiveCompHe3::FiveCompHe3 (int OP_size, complex<double> *etab) : SC_op (OP_size, etab) {};
// --- need to add default constructor ... SC_op () {}; ????
//destructor : automatically calls ~SC_op() ?
FiveCompHe3::~FiveCompHe3() {} 

double FiveCompHe3::bulkRHS(double t, double p, complex<double>* eta, complex<double>* dFdeta, double& FEbulk) {
	Matrix3cd A;
	A<< eta[0],  0 ,  eta[3],
		0,     eta[1],  0,
		eta[4],  0,   eta[2];
	Matrix3cd dFdA;

	he3bulk(t, p, A, dFdA, FEbulk); 

	dFdeta[0] = dFdA(0,0);
	dFdeta[1] = dFdA(1,1);
	dFdeta[2] = dFdA(2,2);
	dFdeta[3] = dFdA(0,2);
	dFdeta[4] = dFdA(2,0);

	return 0.0;
} 
