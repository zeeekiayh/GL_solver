#ifndef _includes
#define _includes

#include <iostream>
#include <fstream> // for file in/out
#include <string>
#include "math.h"
#include <complex>
#include <vector>
#include <chrono> // for timing
#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>
#include <eigen/unsupported/Eigen/KroneckerProduct> // for the Kronecker product
#include "ConvergenceAccelerator.hpp"

using std::cout;
using std::endl;
using std::cin;
using std::complex;
using std::ostream;
using std::vector;
using std::string;
using namespace Eigen;
// using namespace std::complex_literals; // for easy comlpex notation

typedef Triplet<double> Tr;
typedef SparseMatrix<double> SpMat_d;
typedef SparseMatrix<dcomplex> SpMat_cd;
// typedef Triplet<dcomplex> cTr; // will we actually use this?

// returns the unique id corresponding to each op-component in the mesh.
//    n_u: mesh index along u
//    n_v: mesh index along v
//    i:   vector (OP) index
//    size: num elements in mesh
int ID(int size, int n_u, int n_u_max, int n_v, int i) { return size*i + n_u_max*n_v + n_u; }

double abs2(dcomplex x) { return pow(abs(x),2); }
double abs2(double x) { return pow(abs(x),2); }

// // To see a small portion of a (sparse) matrix; for debugging
// void Matrix_SubView(SpMat_d matrix, int n_u, int n_v, int width, int height)
// {
//    for (int h = 0; h < height; h++) {
//       for (int w = 0; w < width; w++) {
//          cout << matrix.coeffRef(n_u+w, n_v+h);
//          if (w+1 < width) cout << ", ";
//       }
//       cout << endl;
//    }
// }

// Insert the matrix "spMat" into the location (i,j) in a
//    sparse matrix of size "size" and return the matrix
template<typename Scalar_type>
SparseMatrix<Scalar_type> Place_subMatrix(int i, int j, int size, SparseMatrix<Scalar_type> spMat)
{
   SparseMatrix<Scalar_type> Mij(size,size);
   Mij.insert(i,j) = 1.;
   return kroneckerProduct(Mij,spMat);
}

// constants used in the GL equations
struct GL_param { double K1, K2, K3, B1, B2, B3, B4, B5, alpha; };

// values used for set up, initial conditions, and algorithmic parameters
struct in_conditions
{
   int SIZEu;    // the number of points...
   int SIZEv;    // "" ...
   int SIZEw;    // "" in orthogonal directions
   int N_loop;   // count limit for iterations
   int maxStore; // max num vectors kept in accerleration
   int wait;     // iterations to wait before using acceleration
   double ACCUR; // the minimum desired accuracy
   double STEP;  // step size
   double rel_p; // relaxation param
};

// ========================================================
// A class for the Order Parameter component; just at one 
//    point on the mesh. Scalar_type must be like dcomplex,
//    double, or similar. If your order parameter is a
//    matrix, it will be flattened as a VectorXcd.
   template <typename Scalar_type>
   class OrderParam
   {
      protected:
      // store it as a vector; the OP at one point
      int num_comp = 1; // we'll assume it's 1, unless in the derived class

      private:
      Matrix<Scalar_type,-1,1> OP;

      public:
      OrderParam() {}
      OrderParam(int n): num_comp(n) {}

      void initialize(int n)
      {
         num_comp = n;
         if (num_comp > 1) OP.resize(num_comp);
      }

      // void Set_OP(Matrix<Scalar_type,-1,1> op) { OP = op; }
      void Set_OP(Matrix<Scalar_type,-1,1> op)
      {
         if (op.rows() != num_comp)
         {
            cout << "ERROR: vector given to Set_OP is not the correct size." << endl;
            return;
         }
         OP = op;
      }

      Scalar_type& operator() (int i)
      {
         // cout << "&operator(i); i = " << i << endl;
         if (i > num_comp-1) // check the range first
            throw "ERROR: index out of range OrderParam::operator()\n";
         // cout << "OP = " << OP << endl;
         return OP(i); // it's component
      }

      OrderParam& operator=(OrderParam& rhs)
      {
         OP = rhs.OP;
         num_comp = rhs.num_comp;

         return *this;
      }
      
      // get the op components into a vector form from a 3x3 matrix
      // ...do we actually need this?
      void Set_OP(Matrix<Scalar_type,3,3> op)
      {
         // flatten the matrix (row major)
         int i = 0; // count the values put into the vector OP

         for (int a = 0; a < 3; a++) {          // For each spin index...
            for (int j = 0; j < 3; j++) {       //   go across all orbital indexes
               if (abs(op(a,j)) > pow(10,-8)) { // If not effectively 0
                  if (i > num_comp) cout << "WARNING: more elements in matrix than specified by num_comp." << endl;
                  else {
                     OP(i) = op(a,j);
                     i++;
                  }
               }
            }
         } // for's
      }

      Matrix<Scalar_type,3,3> GetMatrixForm()
      {
         Matrix<Scalar_type,3,3> mat;
         switch (num_comp)
         {
         case 1:
            // add desired OP forms here
            break;
            // ...
         case 3:
            mat << OP(0), 0.,    0.,
                   0.,    OP(1), 0.,
                   0.,    0.,    OP(2);
            break;
         case 5:
            mat << OP(0), 0.,    OP(1),
                   0.,    OP(2), 0.,
                   OP(3), 0.,    OP(4);
         default:
            cout << "ERROR: Calling function for OP matrix form, but the form for num_comp = " << num_comp << "has not been defined." << endl;
            break;
         }
         return mat;
      }

      dcomplex F_Bulk(GL_param gl)
      {
         // calculate the used forms of A
         auto A = GetMatrixForm();
         auto A_tran = A.transpose();
         auto A_conj = A.conjugate();
         auto A_dag  = A.adjoint();

         double Beta_B = gl.B1+gl.B2 + (gl.B3+gl.B4+gl.B5)/3.; // the bulk beta value

         return -( gl.B1*pow( abs((A * A_tran).trace()), 2)
                  +gl.B2*pow( (A * A_dag).trace(), 2)
                  +gl.B3*(A * A_tran * A_conj * A_dag).trace()
                  +gl.B4*(A * A_dag * A * A_dag).trace()
                  +gl.B5*(A * A_dag * A_conj * A_tran).trace()
                  )/(Beta_B * 9.)
                  +2./3.*(A * A_dag).trace();
         // TODO: return just the real part...
      }
   };
// ========================================================


// ===========================================
// Boundary condition structure for a sigle OP
   struct Bound_Cond
   {
      string typeB, typeT, typeL, typeR; // type of BC: Dirichlet (value) or Neumann (derivative)
      double bB, bT, bL, bR; // for Neumann BC:   slip length
                             // for Dirichlet BC: function value
      Bound_Cond& operator= (Bound_Cond& rhs)
      {
         typeB = rhs.typeB;
         typeT = rhs.typeT;
         typeL = rhs.typeL;
         typeR = rhs.typeR;
         bB = rhs.bB;
         bT = rhs.bT;
         bL = rhs.bL;
         bR = rhs.bR;
         return *this;
      }
   };
// ===========================================


// Access OP elements from a vector that has been flatened from a
//   matrix of size 'sizeU' X 'sizeV', with an OP size 'OPsize'
template<typename Scalar_type>
OrderParam<Scalar_type> matrix_operator(Matrix<Scalar_type,-1,1> vec, int u, int v, int sizeU, int sizeV, int OPsize)
{
   OrderParam<Scalar_type> op(OPsize); // initialize the OP to return
   for (int vi = 0; vi < OPsize; vi++) op(vi) = vec(ID(sizeU*sizeV,u,sizeU,v,u)); // insert all OP components into the OrderParam object
   return op;
}


// ====================================================
// A class that holds the mesh of OP components, builds
//    the whole problem (matrices, rhs vector), and
//    solves the system using (accelerated) relaxation
   template <class Scalar_type>
   class GL_Solver
   {
      // VARIABLES
      private:

      protected:
      int size; // number of mesh points (size of the D matrices or OP-component vectors)
      int OP_size; // number of OP components
      bool update; // says whether or not to use the no_update vector
      GL_param gl; // temperature-dependent parameters for the GL equation
      string method; // if Solve will use normal relaxtion or the accelerated
      in_conditions cond; // struct of all the BC's and other parameters for the methods
      vector<int> no_update; // stores all the indeces that will not be modified in the RHS
      Matrix<Scalar_type,-1,1> solution; // to store the solution to the GL equ. (in the single vector form)
      Matrix<Scalar_type,-1,1> op_vector; // the vector form of the op_matrix
      SparseMatrix<Scalar_type> SolverMatrix; // solver matrix
      SparseMatrix<Scalar_type> Du2, Dw2, Duw; // derivative matrices
      Bound_Cond Auu_BC,Auw_BC,Avv_BC,Awu_BC,Aww_BC; // boundary conditions for OP components
      Matrix<OrderParam<Scalar_type>,-1,-1> op_matrix; // the matrix of OP at each mesh point

      public:
      // CONSTRUCTORS & DECSTRUCTOR
      GL_Solver() {};
      GL_Solver(string conditions_file)
      {
         ReadConditions(conditions_file);
         size = cond.SIZEu*cond.SIZEv;
      }
      // ~GL_Solver() {};

      // virtual functions; to be defined in derived classes.
      //   these 3 must be pure virtual so that they can
      //   be called by other functions in this class
      virtual SparseMatrix<Scalar_type> BuildSolverMatrix() = 0;
      virtual Matrix<Scalar_type,-1,1> makeGuess(Matrix<Scalar_type,-1,1>&) = 0;
      virtual Matrix<Scalar_type,-1,1> RHS() = 0;
      // non-virtual functions; to be defined in derived classes
      Scalar_type F_Grad(int, int, int);
      double free_energy();

      // METHODS
      Matrix<Scalar_type,-1,1> getSolution() const { return solution; }
      int getSolverMatrixSize() const { return SolverMatrix.cols(); }
      in_conditions getConditions() const { return cond; }

      // read in the conditions from the file
      void ReadConditions(string conditions_file)
      {
         string line;
         std::ifstream conditions(conditions_file);

         // get conditions from the file
         if (conditions.is_open()) {
            while (!conditions.eof()) {
               getline(conditions,line);
               if (line[0] == '#') {           // any good way for error handling here?
                  string ls = line.substr(1);
                  if (ls == string("SIZE"))
                  {
                     // the sizes can never be 0, but we'll set them to one if = 0
                     conditions >> cond.SIZEu;
                     if (!cond.SIZEu) cond.SIZEu = 1;
                     conditions >> cond.SIZEv;
                     if (!cond.SIZEv) cond.SIZEv = 1;
                  }
                  else if (ls == string("METHOD"))    conditions >> method;
                  else if (ls == string("update"))    conditions >> update;
                  else if (ls == string("ACCUR"))     conditions >> cond.ACCUR;
                  else if (ls == string("h"))         conditions >> cond.STEP;
                  else if (ls == string("num loops")) conditions >> cond.N_loop;
                  else if (ls == string("rel_p"))     conditions >> cond.rel_p;
                  else if (ls == string("maxStore"))  conditions >> cond.maxStore;
                  else if (ls == string("wait"))      conditions >> cond.wait;
                  else if (ls == string("betas"))
                  {
                     conditions >> gl.B1;
                     conditions >> gl.B2;
                     conditions >> gl.B3;
                     conditions >> gl.B4;
                     conditions >> gl.B5;
                  }
                  else if (ls == string("alpha")) conditions >> gl.alpha;
                  else if (ls == string("Ks"))
                  {
                     conditions >> gl.K1;
                     conditions >> gl.K2;
                     conditions >> gl.K3;
                  }
               }
            }
            conditions.close();
            cout << "NOTICE: using " << method << " to solve." << endl;
         }
         else cout << "Unable to open file:" << conditions_file << endl;
      }

      void ReadBoundaryConditions(string boundary_conditions_file)
      {
         string line;
         std::ifstream BCs(boundary_conditions_file);

         // get boundary conditions from the file
         if (BCs.is_open()) {

            while (line != "#OP size") getline(BCs,line); // find the line with the size
            BCs >> OP_size;
            
            while (!BCs.eof()) {
               getline(BCs,line);
               if (line[0] == '#') {           // any good way for error handling here?
                  string ls = line.substr(1);

                  // Auu
                  if (ls == string("Axx bTop")) {
                     BCs >> Auu_BC.bT;
                     BCs >> Auu_BC.typeT;
                  } else if (ls == string("Axx bBott")) {
                     BCs >> Auu_BC.bB;
                     BCs >> Auu_BC.typeB;
                  }

                  // Auw
                  if (ls == string("Axz bTop")) {
                     BCs >> Auw_BC.bT;
                     BCs >> Auw_BC.typeT;
                  } else if (ls == string("Axz bBott")) {
                     BCs >> Auw_BC.bB;
                     BCs >> Auw_BC.typeB;
                  }

                  // Avv
                  else if (ls == string("Ayy bTop")) {
                     BCs >> Avv_BC.bT;
                     BCs >> Avv_BC.typeT;
                  } else if (ls == string("Ayy bBott")) {
                     BCs >> Avv_BC.bB;
                     BCs >> Avv_BC.typeB;
                  }

                  // Awu
                  if (ls == string("Azx bTop")) {
                     BCs >> Awu_BC.bT;
                     BCs >> Awu_BC.typeT;
                  } else if (ls == string("Azx bBott")) {
                     BCs >> Awu_BC.bB;
                     BCs >> Awu_BC.typeB;
                  }

                  // Aww
                  else if (ls == string("Azz bTop")) {
                     BCs >> Aww_BC.bT;
                     BCs >> Aww_BC.typeT;
                  } else if (ls == string("Azz bBott")) {
                     BCs >> Aww_BC.bB;
                     BCs >> Aww_BC.typeB;
                  }
               }
            }

            BCs.close();
         }
         else cout << "Unable to open file:" << boundary_conditions_file << endl;
      }
   
      // Build the derivative matrices
      void Build_D_Matrices()
      {
         // make vectors to hold the triplets of coefficients
         vector<Tr> coeffs_u2, coeffs_v2, coeffs_uv; // coeffs_y2, coeffs_xy, coeffs_yz

         // u, v, w are orthogonal basis for the system
         //   the n's are indexes in each direction

         // Loop through the entire mesh--row/col order does not matter here
         for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
            for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
               int id = ID(size,n_u,cond.SIZEu,n_v,0);
               
               // InsertCoeff_Du2(id, n_u,   n_v, -2., coeffs_u2);
               // InsertCoeff_Du2(id, n_u-1, n_v,  1., coeffs_u2);
               // InsertCoeff_Du2(id, n_u+1, n_v,  1., coeffs_u2);
               
               InsertCoeff_Dw2(id, n_u, n_v,  -2., coeffs_v2);
               InsertCoeff_Dw2(id, n_u, n_v-1, 1., coeffs_v2);
               InsertCoeff_Dw2(id, n_u, n_v+1, 1., coeffs_v2);
               
               // InsertCoeff_Duw(id, n_u-1, n_v-1, 1./4., coeffs_uv);
               // InsertCoeff_Duw(id, n_u+1, n_v-1,-1./4., coeffs_uv);
               // InsertCoeff_Duw(id, n_u-1, n_v+1,-1./4., coeffs_uv);
               // InsertCoeff_Duw(id, n_u+1, n_v+1, 1./4., coeffs_uv);
               
            }
         }

         // initialize the D's by size
         // Du2.resize(size,size);
         Dw2.resize(size,size);
         // Duw.resize(size,size);

         // build all the D's from their coefficient triplet-vectors
         // Du2.setFromTriplets(coeffs_u2.begin(), coeffs_u2.end());
         Dw2.setFromTriplets(coeffs_v2.begin(), coeffs_v2.end());
         // Duw.setFromTriplets(coeffs_uv.begin(), coeffs_uv.end());
      }

      // // Insert method for the Dx^2 matrix derivatives
      // void InsertCoeff_Du2(int id, int u, int v, double weight, vector<Tr>& coeffs)
      // {
      //    // id is the index of the connecting element, so
      //    //   id1 is the index of the element being inserted
      //    int id1 = ID(size,u,cond.SIZEv,v,0);
      //    // would add boundary conditions here, but we'll use ghost points, so do nothing
      //    if (u == -1 || u == cond.SIZEu);
      //    // we'll just insert default values, filling the spaces,
      //    //   but they will be modified later for each OP-component
      //    else coeffs.push_back(Tr(id,id1,weight));
      // }

      // insert method for the Dz^2 matrix derivatives
      void InsertCoeff_Dw2(int id, int u, int v, double weight, vector<Tr>& coeffs)
      {
         if (v == -1 || v == cond.SIZEv){} // would add boundary conditions here, but
         else coeffs.push_back(Tr(id,ID(size,u,cond.SIZEu,v,0),weight)); // we'll use ghost points, so do nothing
      }

      // // insert method for the mixed derivative matrices
      // void InsertCoeff_Duw(int id, int u, int v, double weight, vector<Tr>& coeffs)
      // {
      //         if (u == -1 || u == cond.SIZEu){} // would add boundary conditions here,
      //    else if (v == -1 || v == cond.SIZEv){} //  but we'll use ghost points, so do nothing
      //    else coeffs.push_back(Tr(id,ID(size,u,cond.SIZEu,v,0),weight));
      // }
   
      SparseMatrix<Scalar_type> Du2_BD(Bound_Cond BC, int op_elem_num)
      {
         SparseMatrix<Scalar_type> Du2_copy;
         return Du2_copy;
      }

      // derivative matrix methods
      SparseMatrix<Scalar_type> Dw2_BD(Bound_Cond BC, int op_elem_num)
      {
         vector<int> indexes_to_visit; // vector for debugging
         SparseMatrix<Scalar_type> Dw2_copy = Dw2;// the matrix that we will edit and return to not modify the original

         int sz = cond.SIZEu*cond.SIZEv; // size of D matrix (num of mesh points)

         for (int n_u = 0; n_u < cond.SIZEu; n_u++) // loop through just the top and bottom boundary points of the mesh
         {
            // indexes for the points on the bottom side
            int id0 =         ID(sz, n_u, cond.SIZEu, 0, 0),
                id0_connect = ID(sz, n_u, cond.SIZEu, 1, 0);
            
            // indexes for the points on the top side
            int idN =         ID(sz, n_u, cond.SIZEu, cond.SIZEv-1, 0),
                idN_connect = ID(sz, n_u, cond.SIZEu, cond.SIZEv-2, 0);

            // set the values at these indexes using the ghost points,
            //   and depending on what kind of BC we have there
            if (BC.typeB == string("Neumann"))
            {
               Dw2_copy.coeffRef(id0,id0) = -2. -2.*cond.STEP/BC.bB;
               Dw2_copy.coeffRef(id0,id0_connect) = 2.;
            }
            else if (BC.typeB == string("Dirichlet"))
            {
               Dw2_copy.coeffRef(id0,id0) = 1.;
               if (!update) no_update.push_back(ID(sz,n_u,cond.SIZEu,0,op_elem_num));
               Dw2_copy.coeffRef(id0,id0_connect) = 0.; // make sure to disconnect from the other connection
            }

            if (BC.typeT == string("Neumann"))
            {
               Dw2_copy.coeffRef(idN,idN) = -2. +2.*cond.STEP/BC.bT;
               Dw2_copy.coeffRef(idN,idN_connect) = 2.;
            }
            else if (BC.typeT == string("Dirichlet"))
            {
               Dw2_copy.coeffRef(idN,idN) = 1.;
               if (!update) no_update.push_back(ID(sz,n_u,cond.SIZEu,cond.SIZEv-1,op_elem_num));
               Dw2_copy.coeffRef(idN,idN_connect) = 0.; // make sure to disconnect from the other connection
            }

            // // debugging (for a large matrix):
            // if (!(n_u%21-2)) indexes_to_visit.push_back(id0);
         }

         // // debugging (for a large matrix):
         // cout << endl << "ouside loop:";
         // for (auto it = indexes_to_visit.begin(); it != indexes_to_visit.end(); it++)
         // {
         //    cout << endl << "mini view:" << endl;
         //    Matrix_SubView(Dw2_copy,*it-2,*it-2,7,7);
         // }

         return Dw2_copy;
      }

      SparseMatrix<Scalar_type> Duw_BD(Bound_Cond BC, int op_elem_num)
      {
         SparseMatrix<Scalar_type> Duw_copy;
         return Duw_copy;
      //    // ?? We actually will not need to add anything to the no_update vector
      //    //   because we have already gone through all the boundary points.
      //    // the matrix that we will edit and return to not modify the original
      //    SparseMatrix<Scalar_type> Duw_copy = Duw;
      //    int sz = cond.SIZEu*cond.SIZEv; // size of D matrix (num of mesh points)
      //    // loop through just the boundary points of the mesh
      //    for (int n_v = 1; n_v < cond.SIZEv-1; n_v++) // loop through the left and right boundary points of the mesh
      //    {
      //       // indexes for the left side
      //       int id0 =             ID(sz, 0, cond.SIZEu, n_v,   0), // the index of the mesh point we're at
      //           id0_connectB =    ID(sz, 0, cond.SIZEu, n_v-1, 0), // index of the bottom point to connect to
      //           id0_connectT =    ID(sz, 0, cond.SIZEu, n_v+1, 0), // index of the top point to connect to
      //          // we need to disconnect from the points that the default D matrix has
      //           id0_disconnectB = ID(sz, 1, cond.SIZEu, n_v-1, 0), // index of the bottom point to disconnect from
      //           id0_disconnectT = ID(sz, 1, cond.SIZEu, n_v+1, 0); // index of the top point to disconnect from
      //       // indexes for the right side
      //       int idN =             ID(sz, cond.SIZEu-1, cond.SIZEu, n_v,   0), // the index of the mesh point we're at
      //           idN_connectB =    ID(sz, cond.SIZEu-1, cond.SIZEu, n_v-1, 0), // index of the bottom point to connect to
      //           idN_connectT =    ID(sz, cond.SIZEu-1, cond.SIZEu, n_v+1, 0), // index of the top point to connect to
      //          // we need to disconnect from the points that the default D matrix has
      //           idN_disconnectB = ID(sz, cond.SIZEu-2, cond.SIZEu, n_v-1, 0), // index of the bottom point to disconnect from
      //           idN_disconnectT = ID(sz, cond.SIZEu-2, cond.SIZEu, n_v+1, 0); // index of the top point to disconnect from
      //       // set the values at these indexes using the ghost points
      //       if (BC.typeL == string("Neumann"))
      //       {
      //          Duw_copy.coeffRef(id0,id0) = 0.; // disconnect from the point itself
      //          Duw_copy.coeffRef(id0,id0_connectT) = cond.STEP/(2.*BC.bL);
      //          // Duw_copy.coeffRef(id0,id0_connectT) = cond.STEP/(2.*BC.bL) (BC.bL >= pow(10,-7) ? cond.STEP/(2.*BC.bL) : 0);
      //          Duw_copy.coeffRef(id0,id0_connectB) = -cond.STEP/(2.*BC.bL);
      //       }
      //       else if (BC.typeL == string("Dirichlet")) Duw_copy.coeffRef(id0,id0) = 1.;
      //       if (BC.typeR == string("Neumann"))
      //       {
      //          Duw_copy.coeffRef(idN,idN) = 0.; // disconnect from the point itself
      //          Duw_copy.coeffRef(idN,idN_connectT) = cond.STEP/(2.*BC.bR);
      //          Duw_copy.coeffRef(idN,idN_connectB) = -cond.STEP/(2.*BC.bR);
      //       }
      //       else if (BC.typeR == string("Dirichlet")) Duw_copy.coeffRef(idN,idN) = 1.;
      //       // disconnect from default connections
      //       Duw_copy.coeffRef(id0,id0_disconnectB) = 0.;
      //       Duw_copy.coeffRef(id0,id0_disconnectT) = 0.;
      //       Duw_copy.coeffRef(idN,idN_disconnectB) = 0.;
      //       Duw_copy.coeffRef(idN,idN_disconnectT) = 0.;
      //       // If we know the VALUE at the point, add the index of it of the guess/solution
      //       //   vector that we already know, i.e. we don't need to update them
      //       if (BC.typeL == string("Dirichlet")) no_update.push_back(ID(sz,0,cond.SIZEu,n_v,op_elem_num));
      //       if (BC.typeR == string("Dirichlet")) no_update.push_back(ID(sz,cond.SIZEu-1,cond.SIZEu,n_v,op_elem_num));
      //    }
      //    for (int n_u = 1; n_u < cond.SIZEu-1; n_u++) // loop through the top and bottom boundary points of the mesh
      //    {
      //       // indexes for the bottom side
      //       int id0 =             ID(sz, n_u,   cond.SIZEu, 0, 0), // the index of the mesh point we're at
      //           id0_connectL =    ID(sz, n_u-1, cond.SIZEu, 0, 0), // index of the left point to connect to
      //           id0_connectR =    ID(sz, n_u+1, cond.SIZEu, 0, 0), // index of the right point to connect to
      //          // we need to disconnect from the points that the default D matrix has
      //           id0_disconnectL = ID(sz, n_u-1, cond.SIZEu, 1, 0), // index of the left point to disconnect from
      //           id0_disconnectR = ID(sz, n_u+1, cond.SIZEu, 1, 0); // index of the right point to disconnect from
      //       // indexes for the top side
      //       int idN =             ID(sz, n_u,   cond.SIZEu, cond.SIZEv-1, 0), // the index of the mesh point we're at
      //           idN_connectL =    ID(sz, n_u-1, cond.SIZEu, cond.SIZEv-1, 0), // index of the left point to connect to
      //           idN_connectR =    ID(sz, n_u+1, cond.SIZEu, cond.SIZEv-1, 0), // index of the right point to connect to
      //          // we need to disconnect from the points that the default D matrix has
      //           idN_disconnectL = ID(sz, n_u-1, cond.SIZEu, cond.SIZEv-2, 0), // index of the left point to disconnect from
      //           idN_disconnectR = ID(sz, n_u+1, cond.SIZEu, cond.SIZEv-2, 0); // index of the right point to disconnect from
      //       // set the values at these indexes using the ghost points
      //       if (BC.typeB == string("Neumann"))
      //       {
      //          Duw_copy.coeffRef(id0,id0) = 0.; // disconnect from the point itself
      //          Duw_copy.coeffRef(id0,id0_connectR) = cond.STEP/(2.*BC.bB);
      //          Duw_copy.coeffRef(id0,id0_connectL) = -cond.STEP/(2.*BC.bB);
      //       }
      //       else if (BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id0,id0) = 1.;
      //       if (BC.typeT == string("Neumann"))
      //       {
      //          Duw_copy.coeffRef(idN,idN) = 0.; // disconnect from the point itself
      //          Duw_copy.coeffRef(idN,idN_connectR) = cond.STEP/(2.*BC.bT);
      //          Duw_copy.coeffRef(idN,idN_connectL) = -cond.STEP/(2.*BC.bT);
      //       }
      //       else if (BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(idN,idN) = 1.;
      //       // disconnect from default connections
      //       Duw_copy.coeffRef(id0,id0_disconnectL) = 0.;
      //       Duw_copy.coeffRef(id0,id0_disconnectR) = 0.;
      //       Duw_copy.coeffRef(idN,idN_disconnectL) = 0.;
      //       Duw_copy.coeffRef(idN,idN_disconnectR) = 0.;
      //       // If we know the VALUE at the point, add the index of it of the guess/solution
      //       //   vector that we already know, i.e. we don't need to update them
      //       if (BC.typeB == string("Dirichlet")) no_update.push_back(ID(sz,n_u,cond.SIZEu,0,op_elem_num));
      //       if (BC.typeT == string("Dirichlet")) no_update.push_back(ID(sz,n_u,cond.SIZEu,cond.SIZEv-1,op_elem_num));
      //    }
      //    // special case for the corners
      //    int id, id_disconnect;
      //    // Top left
      //    id = ID(sz,0,cond.SIZEu,cond.SIZEv-1,0);
      //    id_disconnect = ID(sz,1,cond.SIZEu,cond.SIZEv-2,0);
      //    Duw_copy.coeffRef(id,id_disconnect) = 0.;
      //         if (BC.typeL == string("Neumann") && BC.typeT == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.bL*BC.bT);
      //    else if (BC.typeL == string("Dirichlet") || BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.;
      //          // TODO: determine if this assumption is correct, that
      //          //    if the function value is given for one side, we
      //          //    don't have to worry about the derivative condition.
      //          //    (x4, below)
      //    // Top right
      //    id = ID(sz,cond.SIZEu-1,cond.SIZEu,cond.SIZEv-1,0);
      //    id_disconnect = ID(sz,cond.SIZEu-2,cond.SIZEu,cond.SIZEv-2,0);
      //    Duw_copy.coeffRef(id,id_disconnect) = 0.;
      //         if (BC.typeR == string("Neumann") && BC.typeT == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.bR*BC.bT);
      //    else if (BC.typeR == string("Dirichlet") || BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; // here...
      //    // Bottom left
      //    id = ID(sz,0,cond.SIZEu,0,0);
      //    id_disconnect = ID(sz,1,cond.SIZEu,1,0);
      //    Duw_copy.coeffRef(id,id_disconnect) = 0.;
      //         if (BC.typeL == string("Neumann") && BC.typeB == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.bL*BC.bB);
      //    else if (BC.typeL == string("Dirichlet") || BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; //...here...
      //    // Bottom right
      //    id = ID(sz,cond.SIZEu-1,cond.SIZEu,0,0);
      //    id_disconnect = ID(sz,cond.SIZEu-2,cond.SIZEu,1,0);
      //    Duw_copy.coeffRef(id,id_disconnect) = 0.;
      //         if (BC.typeR == string("Neumann") && BC.typeB == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.bR*BC.bB);
      //    else if (BC.typeR == string("Dirichlet") || BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; //...and here
      //    // If we know the VALUE at the point, add the index of it of the guess/solution
      //    //   vector that we already know, i.e. we don't need to update them
      //    if (BC.typeL == string("Dirichlet") || BC.typeT == string("Dirichlet")) no_update.push_back(ID(sz,0,        cond.SIZEu,cond.SIZEv-1,op_elem_num));
      //    if (BC.typeR == string("Dirichlet") || BC.typeT == string("Dirichlet")) no_update.push_back(ID(sz,size[0]-1,cond.SIZEu,cond.SIZEv-1,op_elem_num));
      //    if (BC.typeL == string("Dirichlet") || BC.typeB == string("Dirichlet")) no_update.push_back(ID(sz,0,        cond.SIZEu,0,        op_elem_num));
      //    if (BC.typeR == string("Dirichlet") || BC.typeB == string("Dirichlet")) no_update.push_back(ID(sz,size[0]-1,cond.SIZEu,0,        op_elem_num));
      //    return Duw_copy;
      }
   
      void BuildProblem()
      {
         Build_D_Matrices();
         initialize_OP_matrix();
         SolverMatrix = BuildSolverMatrix();
      }

      // Convert the OP matrix, at all mesh points, into a vector
      void setVectorForm()
      {
         for (int vi = 0; vi < OP_size; vi++) {
            for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
               for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
                  op_vector(ID(size,n_u,cond.SIZEu,n_v,vi)) = op_matrix(n_v,n_u)(vi);
               }
            }
         }
      }

      void initialize_OP_matrix()
      {
         op_matrix.resize(cond.SIZEv,cond.SIZEu); // initialize matrix
         op_vector.resize(cond.SIZEu*cond.SIZEv*OP_size); // initialize op_vector, for the whole thing (size = num_of_mesh_points * num_OP_components)

         // initialize elements in 'matrix'
         for (int n_u = 0; n_u < cond.SIZEu; n_u++) for (int n_v = 0; n_v < cond.SIZEv; n_v++) op_matrix(n_v,n_u).initialize(OP_size);

         setVectorForm();// make the vector form available
      }

      // given the next guess, update the op at all points on the mesh
      //   and make it available in it's vector form.
      void update_OP_Vector(Matrix<Scalar_type,-1,1>& new_guess)
      {
         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
               Matrix<Scalar_type,-1,1> op(OP_size); // a vector the size of # of OP components
               for (int i = 0; i < OP_size; i++) op(i) = new_guess(ID(size,n_u,cond.SIZEu,n_v,i));
               op_matrix(n_v,n_u).Set_OP(op);
            }
         }
         setVectorForm();
      }

      // use the relaxation method and Anderson Acceleration to solve
      void Solve(Matrix<Scalar_type,-1,1>& guess)
      {
         cout << "solving..." << endl;
         auto start = std::chrono::system_clock::now();
         Matrix<Scalar_type,-1,1> f = makeGuess(guess), df(guess.size()); // initialize vectors

         if (no_update.size()) cout << "using no_update" << endl;

         // use LU decomposition to solve the system
         SparseLU<SparseMatrix<Scalar_type>, COLAMDOrdering<int> > solver;
         solver.analyzePattern(SolverMatrix); // without this, Eigen throws: Eigen::Matrix<int, -1, 1>; ... Assertion `index >= 0 && index < size()' failed.
         solver.factorize(SolverMatrix);

         // check to see if the solver failed
         if (solver.info() == Eigen::Success) cout << "\tSolver: successfully built" << endl;
         else if (solver.info() == Eigen::NumericalIssue) // for debugging non-invertable matrices
         {
            cout << "Solver: numerical issues" << endl;
            // for (int k=0; k < SolverMatrix.outerSize(); ++k) for (SparseMatrix<Scalar_type>::InnerIterator it(SolverMatrix,k); it; ++it) cout << "(" << it.row() << "," << it.col() << ")\t";
            return;
         }

         // time to prepare solving method
         auto end = std::chrono::system_clock::now();
         auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
         cout << "\ttime: " << elapsed.count() << " seconds." << endl << endl;

         update_OP_Vector(f); // prepare the matrix for RHS

         int cts = 0; // count loops
         double err;  // to store current error
         Matrix<Scalar_type,-1,1> rhs = RHS(); // the right hand side

         // the acceleration object
         converg_acceler<Matrix<Scalar_type,-1,1>> Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update);
         
         // loop until f converges or until it's gone too long
         do { // use relaxation

            df = solver.solve(rhs)-f; // find the change in f

            if (method == string("acceleration")) Con_Acc.template next_vector<Matrix<Scalar_type,-1,-1>>(f,df,err); // use Anderson Acceleration to converge faster
            else if (method == string("relaxation")) // use normal relaxation
            {
               f += cond.rel_p*df;
               err = df.norm()/f.norm();
            }
            else { cout << "ERROR: Unknown method type given." << endl; return; }

            update_OP_Vector(f); // update the matrix for RHS
            rhs = RHS();     // update rhs
            cts++;           // increment counter

            // for debugging: to see if the solution is oscillating rather than converging
            // for (int i = 0; i < f.size(); i++) if ((i+1)%cond.SIZEu==0) cout << "\tf(" << i << ") = " << f(i) << endl;

            // output approx. percent completed
            cout << "\033[A\33[2K\r" << "estimated: " << round((cts*100.)/cond.N_loop) << "% done" << endl;

         } while(err > cond.ACCUR && cts < cond.N_loop);

         if (err < cond.ACCUR) cout << "Found solution:" << endl;
         else cout << "Result did not converge satifactorily:" << endl;
         cout << "\titerations = " << cts << endl;
         cout << "\trelative error = " << err << endl;

         solution = f;
      }

      // Write all components of the OP, all into one file, of the form:
      //             __x__|__y__|_Auu_|_Auv_| ...
      //              ... | ... | ... | ... | ...
      // separates the components from solution...storing real and imag parts ==> up to 18
      void WriteToFile(Matrix<Scalar_type,-1,1>& vec, string file_name)
      {
         std::ofstream data (file_name);
         if (data.is_open())
         {
            for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
               for (int n_u = 0; n_u < cond.SIZEu; n_u++) {

                  string line = to_string(n_u*cond.STEP) + string("\t")
                              + to_string(n_v*cond.STEP) + string("\t"); // add the position components

                  for (int vi = 0; vi < OP_size; vi++) {
                     dcomplex element = solution(ID(size,n_u,cond.SIZEu,n_v,vi));
                     line += to_string(element.real())                 // add the real and imaginary
                           + string("\t") + to_string(element.imag()); //   components of the solution vector
                     if (vi+1 < size) line += string("\t");
                  }

                  data << line << endl;
               }
            }
         }
         else cout << "Unable to open file: " << file_name << endl;
      }

      void WriteToFile_single_vector(VectorXd& vec, string file_name)
      {
         std::ofstream data (file_name);
         if (data.is_open()) {
            for (int i = 0; i < vec.size(); i++) {
               string line = to_string(i*cond.STEP) + string("\t") + to_string(vec(i));
               data << line << endl;
            }
         }
         else cout << "Unable to open file: " << file_name << endl;
      }

      void WriteToFile_single_vector(VectorXcd& vec, string file_name)
      {
         std::ofstream data (file_name);
         if (data.is_open()) {
            for (int i = 0; i < vec.size(); i++) {
               string line = to_string(i*cond.STEP) + string("\t") + to_string(vec(i).real()) + string("\t") + to_string(vec(i).imag());
               data << line << endl;
            }
         }
         else cout << "Unable to open file: " << file_name << endl;
      }

   }; // GL_solver class

// derived classes for specific OP forms
   template <class Scalar_type>
   class Three_Component_GL_Solver : public GL_Solver<Scalar_type> {};

   template <class Scalar_type>
   class Five_Component_GL_Solver : public GL_Solver<Scalar_type> {};
// ==================================================

// Now include the real and comples class header files,
//   since the parent class are defined
// #include "real_classes.hpp"
#include "complex_classes.hpp"

#endif