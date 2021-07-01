#ifndef _includes
#define _includes

#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
#include </home/izek/Downloads/eigen-3.3.9/Eigen/Dense>
#include </home/izek/Downloads/eigen-3.3.9/Eigen/Sparse>
#include "/home/izek/Downloads/accelerconverg/AcceleratorConverge.hpp"

using std::cout;
using std::endl;
using std::cin;
using namespace Eigen;

struct in_conditions {
   int SIZE;     // N; the number of points
   double ACCUR; // the minimum desired accuracy
   double STEP;  // the step size
   double B;     // the variable parameter
   int N_loop;   // count limit for iterations
};

VectorXd R(VectorXd&, double);
void Write_to_file(VectorXd&, std::string);
VectorXd Relaxation(double[], VectorXd&);
double Free_Energy(double[], VectorXd&);
double Integrand(double, double, double);

VectorXd R(VectorXd& f, double h)
{
   VectorXd f_n = f;
   for (int i = 1; i < f.size()-1; i++) f_n(i) = (std::pow(f(i),3)-f(i));
   f_n(f.size()-1) = 1;
   f_n(0) = 0;
   return f_n;
}

void Write_to_file(VectorXd& f, std::string file)
{
   std::ofstream data (file);
   if (data.is_open())
   {
      data << f;
      data.close();
   }
   else cout << "Unable to open file: " << file << endl;
}

VectorXd Relaxation(double cond[], VectorXd& guess)
{
   SparseMatrix<double> L; // the tri-diagonal matrix
   VectorXd f; // the initial guess
   VectorXd df; // the small change in f

   int s = cond[0];
   double acc = cond[1];
   double h = cond[2];
   double b = cond[3];
   int n_max = cond[4];
   double p = cond[5]; // relaxation factor

   cout << "size s = " << s << endl;
   cout << "step size h = " << h << endl;
   cout << "variable param b = " << b << endl;
   cout << "relaxation param p = " << p << endl;

   L.resize(s,s);
   f.resize(s);
   df.resize(s);

   double h2 = h*h;

   // fill the elements we know
   for (int i = 1; i < s - 1; i++)
   {
      L.insert(i,i-1) = 1./h2;
      L.insert(i,i)   = -2./h2;
      L.insert(i,i+1) = 1./h2;
   }

   // BC: f'(0) = 1/b f(0)
   L.insert(0,0) = -1./h - 1./b;
   L.insert(0,1) = 1./h;

   // BC: f'_N = 0
   // L.insert(s-1,s-2) = -1/h;
   // L.insert(s-1,s-1) = 1/h;

   // BC: f_N = 1
   L.insert(s-1,s-1) = 1.0;

   // fix the ends of f for the BC's
   f(0) = 0;
   f(s-1) = 1.0;
   
   // use LU decomposition to solve the system
   SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
   solver.analyzePattern(L);
   solver.factorize(L);

   // for (int i = 5; i <= 10; i++){ for (int j = 5; j <= 20; j++) {

   // use Anderson Acceleration method (from Anton's code)
   // converg_acceler<VectorXd> AndrAcclr(5, 12, p);

   // loop until f converges or until it's gone too long
   int cts = 0; // count loops
   double err;
   VectorXd rhs;
   do {
      df = solver.solve(R(f,h)) - f;

      f += p*df;
      err = df.norm()/f.norm();

      // AndrAcclr.next_vector<MatrixXd>(f,df,err);
      cts++;
   } while(err > acc && cts < n_max);

   // cout << "(i, j) = (" << i << ", " << j << ")" << endl;
   cout << "cts = " << cts << endl;
   cout << "err = " << err << endl;

   // }}

   // f /= h*h;
   // f(0) *= h;
   // f(s-1) *= h;

   return f;
}

double Free_Energy(double cond[], VectorXd& f)
{
   int s = cond[0];
   double h = cond[2];
   double b = cond[3];

   // |alpha|^2 / 2beta int_{0}^{infty} {dx [ (1-f^2)^2 + 2xi^2 f'^2 ]}
   
   double I = 0, alpha = -1, beta = 1, xi = 1; // normalize constants for now

   double df, dfp;
   for (int i = 0; i < s-1; i++)
   {
      df = (f(i)+f(i+1))/2;  // approx f(i+1/2)
      dfp = (f(i+1)-f(i))/h; // approx f'(i+1/2)
      I += h*Integrand(df, dfp, xi);
   }

   return I*abs(alpha)/(2*beta);
}

double Integrand(double f, double fp, double xi)
{
   return pow(1-f*f,2) + 2*pow(xi*fp,2);
}

#endif