/* =================================================================================================
	Anton Vorontsov
	Montana State University, Jan 2022 

 We are solving self-consistently system of equations 

 	f = RHS(f)    where f is  N-component vector (real or complex) 
 	
	equations for RHS have to be provided separately, and usually it is a complicated function; 
     	solution for f is searched by iterations, but the acceleration makes educated guesses 
	for next vector f based on a stored Nstored previous vectors f

In the main file:

#include "ConvergenceAccelerator.hpp" 

In the source, declare the vectors and their sizes 
	VectorXd f(N), df(N); // the iteration vectors
	VectorXd rhs(N); // right-hand side of the self-consistency
	vector<int> no_update; // provide list of vector components that need not be updated
				// leave empty if all components of f are updated;
				// e.g. if f(3) need to be kept fixed, no_update.push_back( 3 );

	// Initiate the convergence acceleration: declaration and constructor for the convergence acceleration subroutine:
      converg_acceler<VectorXd> AndrAcclr(Nstored, initIter, rlxp, no_update);
      	// parameters (found by experimenting for fastest convergence for a given problem): 
		// Nstored - how many previous vectors to keep in memory (usually 4-10) 
		// initIter - number of simple relaxation iterations in the beginning (1-10) 
		// rlxp - simple relaxation parameter (0.001 - 0.1)
		// no_update - the vector list of vector components that need no updating 

	// provide initial guess for f, e.g.
	for (int i=0; i<N; i++) f(i)= initial guess;

	do{ // iterations
		f_new_guess = RHS(f); // find f_new_guess for current f
		df = f_new_guess-f;  // deviation from the self-consistency, df=0 when self-consistent solution is found
		df = f - f_new_guess;  // better definition for SC OP equation, for relaxation stability 

            AndrAcclr.next_vector<MatrixXd>( f, df, err ); // smart guess
			// on input it sends just calculated f and df, 
            	// and on output it gives the next f=f_new to try, 
			// and the relative error <double> err=norm(f_new-f_input)/norm(f_new);

	} while( some condition on self-consistency convergence, e.g. err > 1e-3);

For complex vectors replace VectorXd, MatrixXd -> VectorXcd,MatrixXcd
=================================================================================================*/

#ifndef ConvergenceAccelerator_hpp_
#define ConvergenceAccelerator_hpp_

#include <complex>
#include <vector>

#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>

//  ------------------ namespaces -------------------------------------------
using namespace std; 
using namespace Eigen; 

// ----- what we are going to keep in memory on each iteration --------------------
template <typename T>
struct keep_in_memory{
	T vec;
	T dev;
	double norm;
	//-------------------- constructors --------------------------------------- //
	keep_in_memory () {}; // default constructor
	// special constructor to be used in push_back
	keep_in_memory(const T v, const T dv, const double n) {vec=v; dev=dv; norm=n;}; // functional initialization
	//keep_in_memory(const T v, const T dv, const double n) : vec(v) , dev(dv), norm(n) {} // member initialization 
};

// ----- the main (template) class for the acceleration ------------------------------------
template <class T>
class converg_acceler{
	protected:
		int MaxStored=2;  // number of vectors used for next guess
		int waitcount=0, waitXiter=1; // start smart guessing on second iteration 
		int Nstored;
		double relax_param=0.1; // simple relaxation parameter
		double coupl_param_prev;
		vector<keep_in_memory<T>> stored; // default constructor called, storage organized 
		T vec_min, dev_min, vec_prev;
		vector<int> no_update_points; 
		bool details=0;
		FILE *record;
	public:

		//constructor 
		converg_acceler () {Nstored=0;}; // default constructor
		converg_acceler (int maxstor, int wtcnt, double rlxp, const vector<int> no_update); 
		//destructor
		~converg_acceler () {};

		int gotstored (void) {return Nstored;}

		template <class U> void next_vector(T & , T &, double & );
		void print_component(unsigned int n);
};






/* ===================================================================================================
 	the .cpp source for the template class functions has to be included into .hpp file,
 	otherwise problems arise at linkage and compilation
 ================================================================================================= */

#include <iostream>
#include <cstdio>
#include <stdlib.h> // srand(), rand()
#include <time.h> // time
#include <iomanip> // Set precisions
#include <numeric> // iota() function for sorting
#include <cmath>
//#include <complex>
//#include <vector>

//-----------------------------------------------------------------------------------------------------------------------

//constructor for base class
template <class T>
converg_acceler<T>::converg_acceler (int mxstr, int wtX, double rlxp, const vector<int> no_update) 
{
	Nstored=0; 
	if ( MaxStored < mxstr ) MaxStored=mxstr;
	waitXiter=wtX;
	relax_param = rlxp;
	no_update_points = no_update;
	// we also create a seed here for the rand() function
	srand (time(NULL)); 
	if(details) record=fopen("accelereator_errors_record.dat", "w");
	return;
}

//------------------------------------- forward defined functions -------------------------------------------------------
// sorting of vectors based on some parameter "x" 
// simplified from https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
// first find permutations p
template <typename T>
std::vector<int> sort_permutation( const std::vector<keep_in_memory<T>>& x)
{
    std::vector<int> p(x.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](int i, int j){ return ( x[i].norm < x[j].norm ); }); // sorting condition
    return p;
}
// and then apply this permutation to arbitrary vectors
template <typename T>
std::vector<T> apply_permutation( const std::vector<T>& x, const std::vector<int>& p)
{
    std::vector<T> sorted_vec(x.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](int i){ return x[i]; });
    return sorted_vec;
}

// -------- the main acceleration code -----------------------------
template <class T>
template <class U>
void converg_acceler<T>::next_vector (T & v, T & dvz, double & error) 
// on input: v,dvz are the latest pair, dvz is deviation from self-consistency for v 
// on output we return the next guess for v
// algorithm of the acceleration is described in M.Eschrig notes 
{
	int i, go_random=0;
	int dimension = (int) v.size();
	//double norm_dev_min_prev=dev_min.norm();
	double guess_coupl, fix_coupl, coupl_param; 
	T devi(dimension);

	if(MaxStored > dimension+1){ 
		MaxStored = dimension+1;
		cout << " \n !!! Anderson Acceleration: reducing number of stored vectors to avoid singular matrices. !!! \n " << endl;
	}

	vec_prev=v; // we'll remember the incoming vector to later find change error on this update 
	double norm_new=dvz.norm();
	if(details) printf("OP update: new vector with deviation =%.4f\n", norm_new);
	waitcount++;

	if(Nstored >= MaxStored) // we reached max storage capacity; and we need to get rid of something
	{
		if( norm_new < 0.1*stored[0].norm) {
			// new vector is a really good guess, so get rid of last half 
			i=(int) MaxStored/2.0;
			stored.erase( stored.end()-i, stored.end() );
		} else if( norm_new > stored[MaxStored-1].norm){ 
			// we went the wrong way completely on the Jump, so remove some bad last vectors 
			i=(int) MaxStored/2.0;
			stored.erase( stored.end()-i, stored.end() );
		} else{
			// if nothing unusual, just get rid of the last element with the largest norm 
			// (sorting will be done at the end of the subroutine)
			stored.pop_back(); 
		}
	}

	// add the currently supplied (v,dvz,norm_new) to the end of the storage vector if everything is OK
	stored.push_back( (keep_in_memory<T>(v,dvz,norm_new)) );

	Nstored=stored.size(); 

	// print out stored norms
	if(details){ printf("norms(dv):\t"); for(i=0; i<Nstored; i++){ printf("(%.4f)  ", stored[i].norm );} printf("\n"); }

	if(details){ printf("Nstored: %d\n", Nstored); }
	if(details){ printf("waitcount: %d\n", waitcount); }
	if(details){ printf("waitXiter: %d\n", waitXiter); }

	// Now we determine the mapping coefficients by minimizing the distance from dev=zero to the dev-hypersurface 
	if( Nstored>1 &&  waitcount > waitXiter ){
		// T and U are the Vector and Matrix classes, should match! 
		// eg. T=VectorXd, U=MatrixXd; or T=Vector4cd, U=Matrix4cd; etc
		T c(Nstored-1), r(Nstored-1); 
		U M(Nstored-1,Nstored-1);

		// the reference vector does not play role for linear equations, we take the one that we just added 
		int ref_v=Nstored-1;

		for ( i=0; i < Nstored; i++ ){ 
			if(i==ref_v) continue;
			devi = stored[i].dev -stored[ref_v].dev;
			for(int j=0; j < Nstored; j++){
				if(j==ref_v) continue;
				M(i,j) = devi.dot( stored[j].dev-stored[ref_v].dev );
			}
			r(i) = - devi.dot( stored[ref_v].dev );
		}

		// find the coefficients that will minimize the norm 
		c = M.colPivHouseholderQr().solve(r); 
		// find the relative error to estimate if the matrix is singular (this should be real small if M is invertible)
		double relative_error = (M*c - r).norm() / r.norm();

		if( relative_error > 1e-6 ){ 
			if(details) cout << " Andrs Accel. has singular matrix. Det(M) = " << M.determinant() << endl; 
			// degeneracy of vectors: the last one is equal to a previous one,
			// or matrix exploded: we eliminate the last one 
			stored.erase( stored.end()-1, stored.end() );
			// ... and will use a random step later on
			go_random=1;
		}
		else{ 
			if(details){ 
				cout << "c = \n" << c.transpose() << endl;
				double avrgnorm=0; for ( i=0; i<Nstored; i++ ) avrgnorm += stored[i].norm/Nstored;
				fprintf(record, "%d \t%e \t%e \t%e \t%e \t%e \n", waitcount, dev_min.norm(), dvz.norm(), avrgnorm, relative_error, abs(M.determinant()) );
				//fprintf(record, "%d \t%e \t%e \t%e \t%e \t%e \n", waitcount, dev_min.norm(), dvz.norm(), avrgnorm, avrgElmM, log10(detM));
			}

			vec_min=stored[ref_v].vec;
			dev_min=stored[ref_v].dev;
			for(i=0; i < Nstored; i++){
				if(i==ref_v) continue;
				vec_min += c(i) * (stored[i].vec - stored[ref_v].vec); 
				dev_min += c(i) * (stored[i].dev - stored[ref_v].dev); 
			}

			devi = vec_min - v;
			guess_coupl = devi.norm()/dvz.norm();

			fix_coupl = 1.0; //norm_dev_min_prev/norm_new ; // correction if we did overshoot the last time
			if (fix_coupl>1) fix_coupl=1;
			
			coupl_param=guess_coupl*fix_coupl;
			coupl_param_prev = coupl_param; 

			if(details) printf("projected min_dev is =%.4f  Jump with = %.4f=(%.4f x %.4f)\n", dev_min.norm(), coupl_param, guess_coupl, fix_coupl);
			v=vec_min; 
			dvz=dev_min; 
		}
	}
	else{
		// we do simple relaxation 
		// the v and dvz are the ones we got on input
		coupl_param_prev= 0.0;
		coupl_param = relax_param;

		if(details) printf("Relaxation =%.4f\n",coupl_param);
		vec_min=v;
		dev_min=dvz;
	}

	Nstored=stored.size(); 

	// special treatment when things are not as good - do a random step
	if(go_random){ 
		// we should use some better vector to get new trial, take least-norm vector 
		v   = vec_min = stored[0].vec;
		dvz = dev_min = stored[0].dev;

		// ... and take a random relaxation step with coupl_param between  (1...6) * relax_param 
		coupl_param = relax_param*(1.0+5.0*rand()/(RAND_MAX));
		coupl_param_prev= 0.0;
		if(details) printf("Random step with = %.4f  \n", coupl_param);
	}

	// next vector to try
	v += coupl_param * dvz;

	// in case we don't want to update some points we return the incoming component
	for(int j=0; j< (int) no_update_points.size(); j++) {
		i=no_update_points[j];
		v(i) = vec_prev(i);
	}

	devi = v-vec_prev; // difference between the new guess and the previous vector 
	error = devi.norm()/v.norm();
	
	// sorting of vectors:
	// first find permutation that will lead to ascending norm of vectors
	std::vector<int> p = sort_permutation(stored); 
	// then apply this permutation to the storage structure
	stored = apply_permutation(stored, p);

	return;
}

template <class T>
void converg_acceler<T>::print_component(unsigned int n)
{
	if( n < stored.size() )
	{
		cout << "component n=" << n << endl;
		cout << "f= \n " << stored[n].vec << endl;
		cout << "df= \n " << stored[n].dev << endl;
		cout << "norm(df)= " << stored[n].norm << endl;
	}
	else
	{
		cout << n << " : NO such component stored\n";
	}
	return;
}


#endif
