#include "structures.hpp" 

using namespace Eigen;

void betas(double t, double p, double *beta);

// bulk values of He3 free energy and of dFdA 3x3 matrix 
// normalized by the bulk values of the He3-B gap and the free energy 
// at temperature t and pressure p
void he3bulk(double t, double p, double & betaB, Eigen::Matrix<T_scalar,3,3> A, Eigen::Matrix<T_scalar,3,3> & dFdA, double & FE)
{
	double beta[6];
	auto Adag = A.adjoint(); // Hermitian conjugate
	auto At = A.transpose(); // transpose 
	auto Ac = A.conjugate(); // complex conjugate

	betas(t,p, beta);
	// beta_B  = beta12 + beta345/3 // ? DO WE NEED THIS ANYMORE ?
	betaB = beta[1]+beta[2] + (beta[3]+beta[4]+beta[5])/3.0;

	// derivative of FE w.r.to A^* -components 
	dFdA = beta[0]*A; 
	dFdA 	+= 1.0/3.0*beta[1] * (A*At).trace() *Ac;
	dFdA 	+= 1.0/3.0*beta[2] * (A*Adag).trace() *A;
	dFdA 	+= 1.0/3.0*beta[3] * (A * At * Ac);
	dFdA 	+= 1.0/3.0*beta[4] * (A * Adag * A);
	dFdA 	+= 1.0/3.0*beta[5] * (Ac * At * A);

	// bulk free energy, normalized by B-phase value
	T_scalar FElocal 
			 = 2.0/3.0 * beta[0] * (A * Adag).trace();
	FElocal 	+= 1.0/9.0 * beta[1] * (A*At).trace() * (Ac*Adag).trace();
	FElocal 	+= 1.0/9.0 * beta[2] * (A*Adag).trace() * (Ac*At).trace();
	FElocal 	+= 1.0/9.0 * beta[3] * (A * At * Ac * Adag ).trace();
	FElocal 	+= 1.0/9.0 * beta[4] * (A * Adag * A * Adag).trace();
	FElocal 	+= 1.0/9.0 * beta[5] * (A * Adag * Ac * At ).trace();
	// FE = real( FElocal ); // for complex-valued solvers // TODO: is there a better way to do this?
	FE = ( FElocal );        // for real-valued solvers    //       it depends on the typedefs in "structures.hpp"

	return;
}

// beta_i/beta_B parameters 
// beta_B  = beta12 + beta345/3
// include (t,p)-dependent strong coupling effects here 
void betas(double t, double p, double *beta)
{
	beta[0]=-1.0; // alpha-coefficient in GL / |alpha| 
	// strong-coupling corrections will change beta[1-5]
	beta[1]=-1.0; // reference beta, followed by weak-coupling relations 
	beta[2]=-2*beta[1]; 
	beta[3]=-2*beta[1]; 
	beta[4]=-2*beta[1]; 
	beta[5]=+2*beta[1]; 

	// testing stong-coupling corrections:
	//	+ T_over_Tc * Delta_beta
	double T_over_Tc = 0.9; // 90% of Tc
	// test p = 15 bar
	beta[1] += T_over_Tc * 1.0;
	beta[2] += T_over_Tc * 0.975;
	beta[3] += T_over_Tc * 0.86;
	beta[4] += T_over_Tc * 0.675;
	beta[5] += T_over_Tc * 0.95;

	double betaB = beta[1]+beta[2] + (beta[3]+beta[4]+beta[5])/3.0;

	for(int i=1; i<=5; i++) beta[i] /= betaB;
	
	return;
}
