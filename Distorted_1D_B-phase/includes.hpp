#ifndef _includes
#define _includes

#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
#include <complex>
#include </home/izek/Downloads/eigen-3.3.9/Eigen/Dense>
#include </home/izek/Downloads/eigen-3.3.9/Eigen/Sparse>
// #include "/home/izek/Downloads/accelerconverg/AcceleratorConverge.hpp"

using std::cout;
using std::endl;
using std::cin;
using namespace Eigen; // for matricies and vectors
using std::complex;
using std::ostream;
using namespace std::complex_literals; // for easy comlpex notation

// constants used in the GL equations
struct consts { double K1, K2, K3, B1, B2, B3, B4, B5, At; };

// values used for set up, initial conditions, and algorithmic limits
struct in_conditions
{
   int SIZE, N_loop; // N; the number of points, count limit for iterations
   double ACCUR, STEP, B, rel_p; // the minimum desired accuracy, step size, variable parameter, relaxation param
   consts c; // constants--do we even need this since the equ's are normalized?
   bool BC; // boundary conditions at the end: slope or value
};

// a structure to hold both the parallel and perpendicular
//   OP's that are in A for the distorted B-phase
struct deltas
{
   VectorXcd para, perp;

   deltas(VectorXcd pa, VectorXcd pe): para(pa), perp(pe)
   {
      if (pa.size() != pe.size()) cout << "WARNING: vectors in object 'deltas' are not the same shape." << endl;
   }
   deltas(int size): para(size), perp(size) { }

   deltas operator+ (const deltas& d) const { return deltas(d.para+this->para, d.perp+this->perp); }
   deltas operator- (const deltas& d) const { return deltas(d.para-this->para, d.perp-this->perp); }
   deltas operator* (const double& d) const { return deltas(d * this->para, d * this->perp); }
   deltas operator/ (const double& d) const { return deltas(this->para / d, this->perp / d); }
   deltas& operator= (const deltas& d)
   {
      this->para = d.para;
      this->perp = d.perp;
      return *this;
   }
   // deltas& operator() (int i)
   // {
   //    //
   // }
   void operator+= (const deltas& d) { this->para += d.para; this->perp += d.perp; }

   double norm() const { return std::max(para.norm(), perp.norm()); }
   double abs() const { return std::sqrt(this->norm()); }
   int size() const
   {
      if (this->para.size() != this->perp.size())
      {
         cout << "ERROR: size not defined for object 'deltas', vectors are not the same shape." << endl;
         return -1;
      }
      return para.size();
   }
};

ostream& operator<< (ostream& os, const deltas& d)
{
   os << "Parallel:\n" << d.para << "\nPerpendicular:\n" << d.perp;
   return os;
}

// function templates
deltas R(deltas&, in_conditions);
void Write_to_file(VectorXcd&, std::string, std::string);
deltas Relaxation(SparseMatrix<complex<double>>&, deltas&, in_conditions);
complex<double> Free_Energy(in_conditions, deltas&);
double Integrand(deltas, deltas, double);
void SetUp(SparseMatrix<complex<double>>&, deltas&, in_conditions);
Matrix<complex<double>,3,3> M_index(Matrix<VectorXcd,3,3>, int);

// the right hand side of the matrix equation
deltas R(deltas& del, in_conditions cond)
{
   if (del.para.size() != del.perp.size()) cout << "WARNING: del.para and del.perp are not the same shape." << endl;
   
   deltas d = del; // make a type 'deltas' the same size as del

   //for (int i = 1; i < del.para.size()-1; i++)
   for (int i = 0; i < del.para.size(); i++) // loop through the whole vector
   {
      // use the equations of the diff. eq. as obtained from the GL equations of the distorted B-phase
      d.para(i) = -1./5.*(2.*pow(del.para(i),2)+pow(del.perp(i),2))*conj(del.para(i)) + 2./5.*(2*pow(abs(del.para(i)),2)+pow(abs(del.perp(i)),2))*del.para(i) + 2./5.*pow(abs(del.para(i)),2)*del.para(i)-del.para(i);
      d.perp(i) = -1./15.*(2.*pow(del.para(i),2)+pow(del.perp(i),2))*conj(del.perp(i)) + 2./15.*(2*pow(abs(del.para(i)),2)+pow(abs(del.perp(i)),2))*del.perp(i) + 2./15.*pow(abs(del.perp(i)),2)*del.perp(i)-del.perp(i)/3.;
   }

   // keep the ends the same fixed values to satisfy boundary conditions
   //  are these correct like this?
   if(cond.BC)
   {
      d.para(del.para.size()-1) = 1;
      d.perp(del.perp.size()-1) = 1;
   }
   else
   {
      d.para(del.para.size()-1) = 0;
      d.perp(del.perp.size()-1) = 0;
   }
   d.para(0) = 0;
   d.perp(0) = 0;

   return d;
}

// save a complex vector to 2 files: real and imaginary parts
void Write_to_file(VectorXcd& f, std::string file_c, std::string file_r)
{
   // the real part
   std::ofstream data_r (file_r);
   if (data_r.is_open())
   {
      for (int i = 0; i < f.size(); i++) data_r << f(i).real() << "\n";
      data_r.close();
   }
   else cout << "Unable to open file: " << file_r << endl;

   // the imaginary part
   std::ofstream data_c (file_c);
   if (data_c.is_open())
   {
      for (int i = 0; i < f.size(); i++) data_c << f(i).imag() << "\n";
      data_c.close();
   }
   else cout << "Unable to open file: " << file_c << endl;
}

// set up the L matrix with the initial conditions
void SetUp(SparseMatrix<complex<double>>& L_pa, SparseMatrix<complex<double>>& L_pe, deltas& d, in_conditions cond)
{
   // display parameters to easily reproduce results
   cout << "size s = " << cond.SIZE << endl;
   cout << "step size h = " << cond.STEP << endl;
   cout << "variable param b = " << cond.B << endl;
   cout << "relaxation param p = " << cond.rel_p << endl;

   // make sure the matrix is the correct size
   L_pa.resize(cond.SIZE,cond.SIZE);
   L_pe.resize(cond.SIZE,cond.SIZE);

   double h2 = cond.STEP*cond.STEP; // h^2 to reduce calculations in the for-loop

   // fill the elements we know
   for (int i = 1; i < cond.SIZE - 1; i++)
   {
      L_pa.insert(i,i-1) = 1./h2;
      L_pa.insert(i,i)   = -2./h2;
      L_pa.insert(i,i+1) = 1./h2;

      L_pe.insert(i,i-1) = 1./h2;
      L_pe.insert(i,i)   = -2./h2;
      L_pe.insert(i,i+1) = 1./h2;
   }

   // make these boundary conditions available based on the conditions .txt?
   // BC: d'(0) = 0                  // BC: d'(0) = 1/cond.B d(0)
   L_pa.insert(0,0) = -1./cond.STEP;// - 1./cond.B;
   L_pa.insert(0,1) = 1./cond.STEP;
   // L_pa.insert(0,0) = 1.;

   L_pe.insert(0,0) = 1.;
   // L_pe.insert(0,1) = ;

   if(cond.BC)// BC: f_N = 1
   {
      L_pa.insert(cond.SIZE-1,cond.SIZE-1) = 1.0;
      L_pe.insert(cond.SIZE-1,cond.SIZE-1) = 1.0;
   }
   else // BC: d'_N = 0
   {
      L_pa.insert(cond.SIZE-1,cond.SIZE-2) = -1./cond.STEP;
      L_pa.insert(cond.SIZE-1,cond.SIZE-1) = 1./cond.STEP;

      L_pe.insert(cond.SIZE-1,cond.SIZE-2) = -1./cond.STEP;
      L_pe.insert(cond.SIZE-1,cond.SIZE-1) = 1./cond.STEP;
   }
}

// use relaxation to find the solution to the linear equation
deltas Relaxation(SparseMatrix<complex<double>>& L_pa, SparseMatrix<complex<double>>& L_pe, deltas& guess, in_conditions cond)
{
   deltas f = guess, df = guess;
   // use LU decomposition to solve the system

   // setup the solver for the parallel part
   SparseLU<SparseMatrix<complex<double>>, COLAMDOrdering<int> > solver_pa;
   solver_pa.analyzePattern(L_pa);
   solver_pa.factorize(L_pa);
   // check to see if the solver failed
   if (solver_pa.info() == Eigen::Success) cout << "Solver_pa: success" << endl;
   else if (solver_pa.info() == Eigen::NumericalIssue) { cout << "Solver_pa: numerical issues" << endl; }
   else { cout << "Solver_pa: other failure...exiting 'Relaxation'" << endl; }

   // setup the solver for the perpendicular part
   SparseLU<SparseMatrix<complex<double>>, COLAMDOrdering<int> > solver_pe;
   solver_pe.analyzePattern(L_pe);
   solver_pe.factorize(L_pe);
   // check to see if the solver failed
   if (solver_pe.info() == Eigen::Success) cout << "Sovler_pe: success" << endl;
   else if (solver_pe.info() == Eigen::NumericalIssue) { cout << "Sovler_pe: numerical issues" << endl; }
   else { cout << "Sovler_pe: other failure...exiting 'Relaxation'" << endl; }

   // loop until f converges or until it's gone too long
   int cts = 0; // count loops
   double err;  // to store current error
   deltas rhs = R(f,cond); // the right hand side

   do { // use relaxation
      df.para = solver_pa.solve(rhs.para) - f.para; // find the change in f
      df.perp = solver_pe.solve(rhs.perp) - f.perp; //   parallel and perpendicular
      f += df*cond.rel_p;       // update f carfully
      err = df.norm()/f.norm(); // calculate the relative error
      rhs = R(f,cond);          // update rhs
      cts++;
   } while(err > cond.ACCUR && cts < cond.N_loop);

   cout << "cts = " << cts << endl;
   cout << "err = " << err << endl;

   return f; // return the last value for f
}

// calculate the free-energy loss by the integral of the energy density
complex<double> Free_Energy(in_conditions cond, deltas& f)
{
   complex<double> I = 0, df, dfp; // start the integral sum at 0; initialize delta f and delta f'
   // double alpha = -1, beta = 1, xi = 1; // normalize constants for now

   // the order parameter matrix
   Matrix<VectorXcd,3,3> OP;
   OP(0,0) = f.para;
   OP(1,1) = f.para;
   OP(2,2) = f.perp;

   Matrix<complex<double>,3,3> A, AT, A_dag, A_conj;

   complex<double> f_bulk;
   complex<double> f_grad;

   for (int i = 0; i < cond.SIZE-1; i++)
   {
      // calculate the used forms of A
      A = M_index(OP,i);
      AT = A.transpose();
      A_dag = AT.conjugate();
      A_conj = A.conjugate();

      f_bulk = cond.c.B1*(
                  pow( abs((A * AT).trace()), 2)
                  - 2.*pow( (A * A_dag).trace(), 2)
                  - 2.*(A * AT * A_conj * A_dag).trace()
                  - 2.*(A * A_dag * A * A_dag).trace()
                  + 2.*(A * A_dag * A_conj * AT).trace()
               ) + cond.c.At*(A * A_dag).trace();
      //f_grad = ;


      // use a rectangular integral approximation,
      //   centered at the midpoints
      //

      // df = (f(i)+f(i+1))/2.; // approx f(i+1/2)
      // dfp = (f(i+1)-f(i))/(double)cond.STEP; // approx f'(i+1/2)
      // I += Integrand(df, dfp, xi)*cond.STEP;
   }

   return I;
}

// returns a matrix of the values at index i from each vector in m
Matrix<complex<double>,3,3> M_index(Matrix<VectorXcd,3,3> m, int i)
{
   Matrix<complex<double>,3,3> r;
   for(int j = 0; j < 3; j++) for(int k = 0; k < 3; k++) r(j,k) = m(j,k)(i);
   return r;
}

/*double Integrand(double f, double fp, double xi)
{
   // what kind of value does the integral give?
   // return pow(1.-f*f,2) + 2.*pow(xi*fp,2);
}

//*/

#endif