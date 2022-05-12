#include <cmath>
#include <eigen/Eigen/Dense>

using namespace Eigen;

void betas(double t, double p, double *beta);

// bulk values of He3 free energy and of dFdA 3x3 matrix 
// normalized by the bulk values of the He3-B gap and the free energy 
// at temperature t and pressure p
void he3bulk(double t, double p, double & betaB, Matrix3cd A, Matrix3cd & dFdA, double & FE)
{
	double beta[6];
	auto Adag = A.adjoint(); // Hermitian conjugate
	auto At = A.transpose(); // transpose 
	auto Ac = A.conjugate(); // complex conjugate

	betas(t,p, beta);
	// beta_B  = beta12 + beta345/3
	betaB = beta[1]+beta[2] + (beta[3]+beta[4]+beta[5])/3.0;

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
	FE 	+= 1.0/9.0 * beta[3] * real( (A * At * Ac * Adag ).trace()  );
	FE 	+= 1.0/9.0 * beta[4] * real( (A * Adag * A * Adag).trace() );
	FE 	+= 1.0/9.0 * beta[5] * real( (A * Adag * Ac * At ).trace()  );

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
