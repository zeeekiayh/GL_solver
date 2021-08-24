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
typedef Triplet<dcomplex> cTr;
typedef SparseMatrix<dcomplex> SpMat_cd;

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

// returns the unique id corresponding to each op-component in the mesh.
//    n_u: mesh index along u
//    n_v: mesh index along v
//    i:   vector (OP) index
//    size: num elements in mesh
int ID(int size, int n_u, int n_u_max, int n_v, int i) { return size*i + n_u_max*n_v + n_u; }

double abs2(dcomplex x) { return pow(abs(x),2); }
double abs2(double x) { return pow(abs(x),2); }

// To see a small portion of a (sparse) matrix; for debugging
void Matrix_SubView(SpMat_d matrix, int n_u, int n_v, int width, int height)
{
   for (int h = 0; h < height; h++) {
      for (int w = 0; w < width; w++) {
         cout << matrix.coeffRef(n_u+w, n_v+h);
         if (w+1 < width) cout << ", ";
      }
      cout << endl;
   }
}

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
         this->num_comp = n;
         if (this->num_comp > 1) OP.resize(this->num_comp);
      }

      // void Set_OP(Matrix<Scalar_type,-1,1> op) { OP = op; }
      void Set_OP(Matrix<Scalar_type,-1,1> op)
      {
         if (op.rows() != this->num_comp)
         {
            cout << "ERROR: vector given to Set_OP is not the correct size." << endl;
            return;
         }
         this->OP = op;
      }

      Scalar_type& operator() (int i)
      {
         // cout << "&operator(i); i = " << i << endl;
         if (i > this->num_comp-1) // check the range first
            throw "ERROR: index out of range OrderParam::operator()\n";
         // cout << "OP = " << OP << endl;
         return OP(i); // it's component
      }

      OrderParam& operator=(OrderParam& rhs)
      {
         this->OP = rhs.OP;
         this->num_comp = rhs.num_comp;

         return *this;
      }
      
      // get the op components into a vector form from a 3x3 matrix
      void Set_OP(Matrix<Scalar_type,3,3> op)
      {
         // flatten the matrix (row major)
         int i = 0; // count the values put into the vector OP

         for (int a = 0; a < 3; a++) {          // For each spin index...
            for (int j = 0; j < 3; j++) {       //   go across all orbital indexes
               if (abs(op(a,j)) > pow(10,-8)) { // If not effectively 0
                  if (i > this->num_comp) cout << "WARNING: more elements in matrix than specified by this->num_comp." << endl;
                  else {
                     this->OP(i) = op(a,j);
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
         auto A = this->GetMatrixForm_He3Defect();
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

// Matrix<Scalar_type,-1,1>

// This class is specific to the OP form with Auu_BC, Aww_BC, Avv_BC along the main diagonal

   // template<typename Scalar_type>
   // class Three_ComponentOrderParam : public OrderParam<Scalar_type>
   // {
   //    public:
   //    Three_ComponentOrderParam() {}
   //    // Three_ComponentOrderParam(int n) { this->initialize(n); }
   //    // gives the 3x3 form of the op for this special form
   //    Matrix<Scalar_type,3,3> GetMatrixForm_He3Defect() // this function is specific to one OP structure
   //    {
   //       Matrix<Scalar_type,3,3> mat;
   //       mat << OP(0), 0.,    0.,
   //              0.,    OP(1), 0.,
   //              0.,    0.,    OP(2);
   //       return mat;
   //    }
   // };

   // template<typename Scalar_type>
   // class Five_ComponentOrderParam : public OrderParam<Scalar_type>
   // {
   //    public:
   //    Five_ComponentOrderParam() {}
   //    // Five_ComponentOrderParam(int n) { this->initialize(n); }
   //    // gives the 3x3 form of the op for this special form
   //    Matrix<Scalar_type,3,3> GetMatrixForm_He3Defect() // this function is specific to one OP structure
   //    {
   //       Matrix<Scalar_type,3,3> mat;
   //       mat << OP(0), 0.,    OP(1),
   //              0.,    OP(2), 0.,
   //              OP(3), 0.,    OP(4);
   //       return mat;
   //    }
   // };
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
         this->typeB = rhs.typeB;
         this->typeT = rhs.typeT;
         this->typeL = rhs.typeL;
         this->typeR = rhs.typeR;
         this->bB = rhs.bB;
         this->bT = rhs.bT;
         this->bL = rhs.bL;
         this->bR = rhs.bR;
         return *this;
      }
   };
// ===========================================


// ====================================================
// A class that holds the mesh of OP components, builds
//    the whole problem (matrices, rhs vector), and
//    solves the system using (accelerated) relaxation
   template <class Scalar_type>
   class GL_Solver
   {
      // VARIABLES
      private:
      // Matrix<OrderParam<Container_type,Scalar_type>,-1,-1> op_matrix; // the matrix of OP at each mesh point
      // Matrix<Scalar_type,-1,1> op_vector; // the vector form of the op_matrix

      protected:
      int size; // number of mesh points (size of the D matrices or OP-component vectors)
      int OP_size; // number of OP components
      bool update; // says whether or not to use the no_update vector
      GL_param gl; // temperature-dependent parameters for the GL equation
      string method; // if Solve will use normal relaxtion or the accelerated
      in_conditions cond; // struct of all the BC's and other parameters for the methods
      vector<int> no_update; // stores all the indeces that will not be modified in the RHS
      Bound_Cond Auu_BC,Auw_BC,Avv_BC,Awu_BC,Aww_BC; // boundary conditions for OP components
      Matrix<Scalar_type,-1,1> solution; // to store the solution to the GL equ. (in the single vector form)
      Matrix<Scalar_type,-1,1> op_vector; // the vector form of the op_matrix
      SparseMatrix<Scalar_type> SolverMatrix; // solver matrix
      SparseMatrix<Scalar_type> Du2, Dw2, Duw; // derivative matrices
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

      // METHODS

      // functions to be defined in specialized derived classes
      void setVectorForm();
      void Solve(Matrix<Scalar_type,-1,1>&);
      void initialize_OP_matrix();
      Matrix<Scalar_type,-1,1> makeGuess(Matrix<Scalar_type,-1,1>&);
      void ReadBoundaryConditions(string);
      void BuildProblem(int,Bound_Cond,Bound_Cond);
      void WriteToFile(string); // Write all components of the OP, all into one file

      Matrix<Scalar_type,-1,1> getSolution() const { return solution; }
      int getSolverMatrixSize() const { return SolverMatrix.cols(); }
      in_conditions getConditions() const { return cond; }

      // Scalar_type f_bulk();
      Scalar_type F_Grad(int, int, int);
      double free_energy();

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
   
   }; // GL_solver class

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