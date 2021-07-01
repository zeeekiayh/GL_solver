#include "includes.hpp"

int main();

VectorXd R(VectorXd& f)
{
   VectorXd f_n = f;
   for (int i = 1; i < f.size()-1; i++) f_n(i) = std::pow(f(i),3)-f(i);
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
   SparseMatrix<double> L;
   VectorXd f1; // the initial guess
   VectorXd df; // the small change in f

   int s = cond[0];
   double acc = cond[1];
   double h = cond[2];
   double b = cond[3];
   int n_max = cond[4];
   double p = cond[5]; // relaxation factor

   L.resize(s,s);
   f1.resize(s);
   df.resize(s);

   double h2 = h*h;

   // fill the elements we know
   for (int i = 1; i < s - 1; i++)
   {
      L.insert(i,i-1) = 1/h2;
      L.insert(i,i)   = -2/h2;
      L.insert(i,i+1) = 1/h2;
   }

   // BC: f'(0) = 1/b f(0)
   L.insert(0,0) = -1/h - 1/b;
   L.insert(0,1) = 1/h;

   // BC: f'_N = 0
   // L.insert(s-1,s-2) = -1/h;
   // L.insert(s-1,s-1) = 1/h;

   // BC: f_N = 1
   L.insert(s-1,s-1) = 1.0;

   // fix the ends of f for the BC's
   f1(0) = 0;
   f1(s-1) = 1.0;
   
   // use LU decomposition to solve the system
   SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
   solver.analyzePattern(L);
   solver.factorize(L);

   // loop until f converges or until it's gone too long
   int cts = 0;    // count loops
   do {
      f1 = solver.solve(f1);
      cts++;
   } while(df.norm()/f1.norm() > acc && cts < n_max);

   cout << "cts = " << cts << endl;
   cout << "acc = " << df.norm()/f1.norm() << endl;

   return f1;
}

double Free_Energy()
{
   //

   return 0;
}