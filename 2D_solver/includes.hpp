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

// Insert the matrix "spMat" into the "location" (i,j) in
//   a sparse matrix of size "size" and return the matrix
template<typename Scalar_type>
SparseMatrix<Scalar_type> Place_subMatrix(int i, int j, int size, SparseMatrix<Scalar_type> spMat) {
   SparseMatrix<Scalar_type> Mij(size,size);
   Mij.insert(i,j) = 1.;
   return kroneckerProduct(Mij,spMat);
}

// constants used in the GL equations
struct GL_param { double K1, K2, K3, B1, B2, B3, B4, B5, alpha, P, T; };

// ===========================================
// Values used for set up, initial conditions,
//   and algorithmic parameters
   struct in_conditions {
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
// ===========================================


// ===========================================
// Boundary condition structure for a sigle OP
   struct Bound_Cond {
      string typeB, typeT, typeL, typeR; // type of BC: Dirichlet (value) or Neumann (derivative)
      double valB, valT, valL, valR; // for Neumann BC:   slip length
                                     // for Dirichlet BC: function value
      Bound_Cond& operator= (Bound_Cond& rhs)
      {
         typeB = rhs.typeB;
         typeT = rhs.typeT;
         typeL = rhs.typeL;
         typeR = rhs.typeR;
         valB = rhs.valB;
         valT = rhs.valT;
         valL = rhs.valL;
         valR = rhs.valR;
         return *this;
      }
   };

   ifstream& operator>> (ifstream& stream, Bound_Cond& BC) {
      stream >> BC.typeT;
      stream >> BC.valT;
      stream >> BC.typeB;
      stream >> BC.valB;
      stream >> BC.typeL;
      stream >> BC.valL;
      stream >> BC.typeR;
      stream >> BC.valR;
      return stream;
   }

   // while (fscanf(ptr,"%*s %*s %s ",buf)==1)
   //      printf("%s\n", buf);

   ostream& operator<< (ostream& out, Bound_Cond& BC) {
      out << "BC.typeT: " << BC.typeT << endl;
      out << "BC.valT:  " << BC.valT << endl;
      out << "BC.typeB: " << BC.typeB << endl;
      out << "BC.valB:  " << BC.valB << endl;
      out << "BC.typeL: " << BC.typeL << endl;
      out << "BC.valL:  " << BC.valL << endl;
      out << "BC.typeR: " << BC.typeR << endl;
      out << "BC.valR:  " << BC.valR << endl;
      return out;
   }
// ===========================================

// ========================================================
// A class for the Order Parameter component; just at one 
//    point on the mesh. Scalar_type must be like dcomplex,
//    double, or similar. If your order parameter is a
//    matrix, it will be flattened as a VectorXcd.
   template <typename Scalar_type>
   class OrderParam {
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
         OP.resize(num_comp);
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

      OrderParam& operator=(OrderParam rhs)
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
         auto Eta_ = GetMatrixForm();
         auto Eta_tran = Eta_.transpose();
         auto Eta_conj = Eta_.conjugate();
         auto Eta_dag  = Eta_.adjoint();

         double Beta_B = gl.B1+gl.B2 + (gl.B3+gl.B4+gl.B5)/3.; // the bulk beta value

         return -( gl.B1*pow( abs((Eta_ * Eta_tran).trace()), 2)
                  +gl.B2*pow( (Eta_ * Eta_dag).trace(), 2)
                  +gl.B3*(Eta_ * Eta_tran * Eta_conj * Eta_dag).trace()
                  +gl.B4*(Eta_ * Eta_dag * Eta_ * Eta_dag).trace()
                  +gl.B5*(Eta_ * Eta_dag * Eta_conj * Eta_tran).trace()
                  )/(Beta_B * 9.)
                  +2./3.*(Eta_ * Eta_dag).trace();
         // TODO: return just the real part...
      }
   };
// ========================================================


// A special function to access OP elements from a vector that
//   has been flattened from a matrix of size 'sizeU' X 'sizeV'
template<typename Scalar_type>
OrderParam<Scalar_type> matrix_operator(Matrix<Scalar_type,-1,1>& vec, int v, int u, int sizeV, int sizeU, int OPsize, bool debug = false) {
   OrderParam<Scalar_type> op;
   op.initialize(OPsize); // initialize the OP to return
   if (debug) cout << "size = " << sizeV*sizeU << endl;
   for (int vi = 0; vi < OPsize; vi++) {
      if (debug) cout << "ID = " << ID(sizeU*sizeV,u,sizeU,v,vi) << endl;
      op(vi) = vec(ID(sizeU*sizeV,u,sizeU,v,vi)); // insert all OP components into the OrderParam object
   }
   return op;
}


// Functions for the temperature- and pressure-dependent GL parameters
double Alpha(double P, double T) {
   // ...
   return 1.0;
}
vector<double> Betas(double P, double T) {
   // ...
   vector<double> b;
   b.push_back(1);
   b.push_back(-2);
   b.push_back(-2);
   b.push_back(-2);
   b.push_back(2);
   return b;
}


// ====================================================
// A class that holds the mesh of OP components, builds
//    the whole problem (matrices, rhs vector), and
//    solves the system using (accelerated) relaxation
   template <class Scalar_type>
   class GL_Solver {
      protected: // for variables and objects that will be used in derived classes

      int size; // number of mesh points (size of the D matrices or OP-component vectors)
      int OP_size; // number of OP components
      GL_param gl; // temperature-dependent parameters for the GL equation
      bool update = 1; // says whether or not to use the no_update vector
      in_conditions cond; // struct of all the BC's and other parameters for the methods
      vector<int> no_update; // stores all the indeces that will not be modified in the RHS
      string method = "acceleration"; // if Solve will use normal relaxtion or the accelerated
      Matrix<Scalar_type,-1,1> solution; // to store the solution to the GL equ. (in the single vector form)
      Matrix<Scalar_type,-1,1> op_vector; // the vector form of the op_matrix
      SparseMatrix<Scalar_type> SolverMatrix; // solver matrix
      SparseMatrix<Scalar_type> Du2, Dw2, Duw; // derivative matrices: ADD D-MATRIX NAMES HERE
      Bound_Cond Eta_uu_BC, Eta_uw_BC, Eta_vv_BC, Eta_wu_BC, Eta_ww_BC; // boundary conditions for OP components

      public:
      // CONSTRUCTORS
      GL_Solver() {};
      GL_Solver(string conditions_file) {
         ReadConditions(conditions_file);
         size = cond.SIZEu*cond.SIZEv;
      }

      // Virtual Functions; to be defined in derived classes.
      //   these 3 must be pure virtual so that they can
      //   be called by other functions in this class
      virtual SparseMatrix<Scalar_type> BuildSolverMatrix() = 0;
      virtual Matrix<Scalar_type,-1,1> makeGuess(Matrix<Scalar_type,-1,1>&) = 0;
      virtual Matrix<Scalar_type,-1,1> RHS() = 0;

      // Non-Virtual Functions; to be defined in derived classes
      Scalar_type F_Grad(int, int, int);
      double free_energy();

      // METHODS
      Matrix<Scalar_type,-1,1> getSolution() const { return solution; }
      int getSolverMatrixSize() const { return SolverMatrix.cols(); }
      in_conditions getConditions() const { return cond; }

      // read in the conditions from the file
      void ReadConditions(string conditions_file) {
         // cout << "reading conditions..." << endl;
         string line;
         char ln [20];
         std::ifstream conditions(conditions_file);

         // get conditions from the file
         if (conditions.is_open()) {
            // while (line != "BOUNDARY CONDITIONS") {getline(conditions,line);/*cout<<"line="<<line<<endl;*/}
            getline(conditions,line); getline(conditions,line); //cout << "line = " << line << endl;
            conditions >> Eta_uu_BC; conditions.ignore(256,'\n');
            conditions >> Eta_vv_BC; conditions.ignore(256,'\n');
            conditions >> Eta_ww_BC; conditions.ignore(256,'\n');
            conditions >> Eta_uw_BC; conditions.ignore(256,'\n');
            conditions >> Eta_wu_BC; conditions.ignore(256,'\n');

            // cout << "boundary conditions read in:" << endl;
            // cout << "Eta_uu_BC\n" << Eta_uu_BC << endl;
            // cout << "Eta_vv_BC\n" << Eta_vv_BC << endl;
            // cout << "Eta_ww_BC\n" << Eta_ww_BC << endl;
            // cout << "Eta_uw_BC\n" << Eta_uw_BC << endl;
            // cout << "Eta_wu_BC\n" << Eta_wu_BC << endl;
            
            conditions >> cond.SIZEu >> cond.SIZEv; conditions.ignore(256,'\n');
            if (!cond.SIZEu) cond.SIZEu = 1; // the sizes can never be 0; set them to 1 if = 0
            if (!cond.SIZEv) cond.SIZEv = 1;

            // if(getline(conditions,ln)) sscanf(ln,"%le",&cond.STEP);
            // cout << "step = " << cond.STEP << endl;

            conditions >> cond.STEP;   conditions.ignore(256,'\n');
            conditions >> OP_size;     conditions.ignore(256,'\n');
            conditions >> gl.T;        conditions.ignore(256,'\n');
            conditions >> gl.P;        conditions.ignore(256,'\n');
            conditions >> cond.ACCUR;  conditions.ignore(256,'\n');
            conditions >> cond.N_loop; conditions.ignore(256,'\n');

            // cout << "size = " << cond.SIZEu << ", " << cond.SIZEv << endl;
            // cout << "step = " << cond.STEP << endl;
            // cout << "OP_size = " << OP_size << endl;
            // cout << "T = " << gl.T << endl;
            // cout << "P = " << gl.P << endl;
            // cout << "initial conditions read in." << endl;

            // read the K matrix
            Matrix<Matrix2d,-1,-1> K; // define the k-matrix: matrix of matrices
            K.resize(OP_size,OP_size);

            Matrix2d K0 = MatrixXd::Zero(2,2); // the basic 0-matrix
            Matrix2d K_temp(2,2); // the temporary k-matrix for importing

            // find the line where the K matrix starts
            while (line != "K MATRIX") {getline(conditions,line);/*cout<<"line="<<line<<endl;*/}

            // get K matrix elements from the file
            for (int row = 0; row < OP_size && !conditions.eof(); row++) {
               getline(conditions,line);

               bool started = false; // to check if we're in a set of parentheses
               int k_index = 0;      // index of the inner k-matrices
               int col = 0;          // column number

               for (int i = 0; i < line.length(); i++) {
                  if (line[i] == '(') started = true; // toggle the started flag when "("
                  if (line[i] == ')') {               // or ")" is reached
                     started = false;
                     k_index = 0; // and reset the index

                     // if it was an empty set of parentheses...
                     if (line[i-1] == '(') K_temp = K0; // add the zero-matrix

                     K(row,col) = K_temp; // if we just got to the end of a set, add the inner K matrix
                     col++;               // increment for the next one's position
                  }
                  if (started) {
                     if (isdigit(line[i])) {         // only if it's a number
                        int r = floor(k_index/2);    // for some reason, we have to define the index
                        int c = k_index%2;           //   values separately, then pass them in!
                        K_temp(r,c) = line[i] - '0'; // add the number to the matrix
                        k_index++;                   // increment for the next location in the matrix
                     }
                  }
               } // end for (row)
            } // end file open
            // cout << "K matrix read in." << endl;

            conditions.close();
            cout << "NOTICE: using " << method << " to solve." << endl;
         }
         else cout << "Unable to open file: " << conditions_file << endl;
      }
   
      // ADD CODE IN THIS FUNCTION TO BUILD ADDITIONAL D-MATRICES
      // Build the derivative matrices
      void Build_D_Matrices() {
         // cout << "build matrices..." << endl;
         // make vectors to hold the triplets of coefficients
         vector<Tr> coeffs_u2, coeffs_w2, coeffs_uw; // coeffs_v2, coeffs_uv, coeffs_vw

         // u, v, w are orthogonal basis for the system
         //   the n's are indexes in each direction

         // Loop through the entire mesh--row/col order does not matter here
         for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
            for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
               int id = ID(size,n_u,cond.SIZEu,n_v,0);
               
               InsertCoeff_Du2(id, n_u,   n_v, -2., coeffs_u2);
               InsertCoeff_Du2(id, n_u-1, n_v,  1., coeffs_u2);
               InsertCoeff_Du2(id, n_u+1, n_v,  1., coeffs_u2);
               
               InsertCoeff_Dw2(id, n_u, n_v,  -2., coeffs_w2);
               InsertCoeff_Dw2(id, n_u, n_v-1, 1., coeffs_w2);
               InsertCoeff_Dw2(id, n_u, n_v+1, 1., coeffs_w2);
               
               InsertCoeff_Duw(id, n_u-1, n_v-1, 1./4., coeffs_uw);
               InsertCoeff_Duw(id, n_u+1, n_v-1,-1./4., coeffs_uw);
               InsertCoeff_Duw(id, n_u-1, n_v+1,-1./4., coeffs_uw);
               InsertCoeff_Duw(id, n_u+1, n_v+1, 1./4., coeffs_uw);
            }
         }

         // cout << "looped successfully" << endl;

         // initialize the D's by size
         Du2.resize(size,size);
         Dw2.resize(size,size);
         Duw.resize(size,size);

         // build all the D's from their coefficient triplet-vectors
         Du2.setFromTriplets(coeffs_u2.begin(), coeffs_u2.end());
         Dw2.setFromTriplets(coeffs_w2.begin(), coeffs_w2.end());
         Duw.setFromTriplets(coeffs_uw.begin(), coeffs_uw.end());

         // cout << "build matrices: done" << endl;
      }

      // Insert method for the Du^2 matrix derivatives
      void InsertCoeff_Du2(int id, int u, int v, double weight, vector<Tr>& coeffs) {
         if (u == -1 || u == cond.SIZEu);
         else coeffs.push_back(Tr(id,ID(size,u,cond.SIZEv,v,0),weight));
      }

      // insert method for the Dw^2 matrix derivatives
      void InsertCoeff_Dw2(int id, int u, int v, double weight, vector<Tr>& coeffs) {
         if (v == -1 || v == cond.SIZEv); // would add boundary conditions here, but
         else coeffs.push_back(Tr(id,ID(size,u,cond.SIZEu,v,0),weight)); // we'll use ghost points, so do nothing
      }

      // insert method for the mixed derivative matrices
      void InsertCoeff_Duw(int id, int u, int v, double weight, vector<Tr>& coeffs) {
              if (u == -1 || u == cond.SIZEu); // would add boundary conditions here,
         else if (v == -1 || v == cond.SIZEv); //  but we'll use ghost points, so do nothing
         else coeffs.push_back(Tr(id,ID(size,u,cond.SIZEu,v,0),weight));
      }
   
      // WRITE ADDITIONAL 'InsertCoeff_D**' METHODS HERE

      // ?still needs fixed? I think it's good
      // derivative matrix: 2nd-order of 1st coordinate (i.e. x)
      SparseMatrix<Scalar_type> Du2_BD(Bound_Cond BC, int op_elem_num) {
         // cout << "\t\t\tDu2_BD()" << endl;
         // vector<int> indexes_to_visit; // vector for debugging
         SparseMatrix<Scalar_type> Du2_copy = Du2;// the matrix that we will edit and return to not modify the original

         for (int n_v = 0; n_v < cond.SIZEv; n_v++) // loop through just the top and bottom boundary points of the mesh
         {
            // indexes for the points on the bottom side
            int id0 =         ID(size, 0, cond.SIZEu, n_v, 0),
                id0_connect = ID(size, 1, cond.SIZEu, n_v, 0);
            
            // indexes for the points on the top side
            int idN =         ID(size, cond.SIZEu-1, cond.SIZEu, n_v, 0),
                idN_connect = ID(size, cond.SIZEv-2, cond.SIZEu, n_v, 0);

            // set the values at these indexes using the ghost points,
            //   and depending on what kind of BC we have there
            if (BC.typeB == string("Neumann"))
            {
               Du2_copy.coeffRef(id0,id0) = -2. -2.*cond.STEP/BC.valB;
               Du2_copy.coeffRef(id0,id0_connect) = 2.;
            }
            else if (BC.typeB == string("Dirichlet"))
            {
               Du2_copy.coeffRef(id0,id0) = 1.;
               if (!update) no_update.push_back(ID(size,0,cond.SIZEu,n_v,op_elem_num));
               Du2_copy.coeffRef(id0,id0_connect) = 0.; // make sure to disconnect from the other connection
            }

            if (BC.typeT == string("Neumann"))
            {
               Du2_copy.coeffRef(idN,idN) = -2. +2.*cond.STEP/BC.valT;
               Du2_copy.coeffRef(idN,idN_connect) = 2.;
            }
            else if (BC.typeT == string("Dirichlet"))
            {
               Du2_copy.coeffRef(idN,idN) = 1.;
               if (!update) no_update.push_back(ID(size,cond.SIZEu-1,cond.SIZEu,n_v,op_elem_num));
               Du2_copy.coeffRef(idN,idN_connect) = 0.; // make sure to disconnect from the other connection
            }

            // // debugging (for a large matrix):
            // if (!(n_v%21-2)) indexes_to_visit.push_back(id0);
         }

         // // debugging (for a large matrix):
         // cout << endl << "ouside loop:";
         // for (auto it = indexes_to_visit.begin(); it != indexes_to_visit.end(); it++)
         // {
         //    cout << endl << "mini view:" << endl;
         //    Matrix_SubView(Du2_copy,*it-2,*it-2,7,7);
         // }

         // cout << "\t\t\tDu2_BD(): success" << endl;
         return Du2_copy;
      }

      // derivative matrix: 2nd-order of 3rd coordinate (i.e. z)
      SparseMatrix<Scalar_type> Dw2_BD(Bound_Cond BC, int op_elem_num) {
         // cout << "\t\t\tDw2_BD" << endl;
         // vector<int> indexes_to_visit; // vector for debugging
         SparseMatrix<Scalar_type> Dw2_copy = Dw2;// the matrix that we will edit and return to not modify the original

         for (int n_u = 0; n_u < cond.SIZEu; n_u++) // loop through just the top and bottom boundary points of the mesh
         {
            // indexes for the points on the bottom side
            int id0 =         ID(size, n_u, cond.SIZEu, 0, 0),
                id0_connect = ID(size, n_u, cond.SIZEu, 1, 0);
            
            // indexes for the points on the top side
            int idN =         ID(size, n_u, cond.SIZEu, cond.SIZEv-1, 0),
                idN_connect = ID(size, n_u, cond.SIZEu, cond.SIZEv-2, 0);

            // set the values at these indexes using the ghost points,
            //   and depending on what kind of BC we have there
            if (BC.typeB == string("Neumann"))
            {
               Dw2_copy.coeffRef(id0,id0) = -2. -2.*cond.STEP/BC.valB;
               Dw2_copy.coeffRef(id0,id0_connect) = 2.;
            }
            else if (BC.typeB == string("Dirichlet"))
            {
               Dw2_copy.coeffRef(id0,id0) = 1.;
               if (!update) no_update.push_back(ID(size,n_u,cond.SIZEu,0,op_elem_num));
               Dw2_copy.coeffRef(id0,id0_connect) = 0.; // make sure to disconnect from the other connection
            }

            if (BC.typeT == string("Neumann"))
            {
               Dw2_copy.coeffRef(idN,idN) = -2. +2.*cond.STEP/BC.valB;
               Dw2_copy.coeffRef(idN,idN_connect) = 2.;
            }
            else if (BC.typeT == string("Dirichlet"))
            {
               Dw2_copy.coeffRef(idN,idN) = 1.;
               if (!update) no_update.push_back(ID(size,n_u,cond.SIZEu,cond.SIZEv-1,op_elem_num));
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

         // cout << "\t\t\tDw2_BD: success" << endl;
         return Dw2_copy;
      }

      // still needs fixed?
      // mixed derivative: of 1st and 3rd coordinates (i.e. x & z)
      SparseMatrix<Scalar_type> Duw_BD(Bound_Cond BC, int op_elem_num) {
         // cout << "\t\t\tDuw_BD" << endl;
         // the matrix that we will edit and return to not modify the original
         SparseMatrix<Scalar_type> Duw_copy = Duw;

         // loop through just the boundary points of the mesh
         for (int n_v = 1; n_v < cond.SIZEv-1; n_v++) // loop through the left and right boundary points of the mesh
         {
            // indexes for the left side
            int id0 =             ID(size, 0, cond.SIZEu, n_v,   0), // the index of the mesh point we're at
                id0_connectB =    ID(size, 0, cond.SIZEu, n_v-1, 0), // index of the bottom point to connect to
                id0_connectT =    ID(size, 0, cond.SIZEu, n_v+1, 0), // index of the top point to connect to
                // we need to disconnect from the points that the default D matrix has
                id0_disconnectB = ID(size, 1, cond.SIZEu, n_v-1, 0), // index of the bottom point to disconnect from
                id0_disconnectT = ID(size, 1, cond.SIZEu, n_v+1, 0); // index of the top point to disconnect from
            
            // indexes for the right side
            int idN =             ID(size, cond.SIZEu-1, cond.SIZEu, n_v,   0), // the index of the mesh point we're at
                idN_connectB =    ID(size, cond.SIZEu-1, cond.SIZEu, n_v-1, 0), // index of the bottom point to connect to
                idN_connectT =    ID(size, cond.SIZEu-1, cond.SIZEu, n_v+1, 0), // index of the top point to connect to
                // we need to disconnect from the points that the default D matrix has
                idN_disconnectB = ID(size, cond.SIZEu-2, cond.SIZEu, n_v-1, 0), // index of the bottom point to disconnect from
                idN_disconnectT = ID(size, cond.SIZEu-2, cond.SIZEu, n_v+1, 0); // index of the top point to disconnect from
            
            // set the values at these indexes using the ghost points
            if (BC.typeL == string("Neumann"))
            {
               Duw_copy.coeffRef(id0,id0) = 0.; // disconnect from the point itself
               Duw_copy.coeffRef(id0,id0_connectT) = cond.STEP/(2.*BC.valL);
               // Duw_copy.coeffRef(id0,id0_connectT) = cond.STEP/(2.*BC.valL) (BC.valL >= pow(10,-7) ? cond.STEP/(2.*BC.valL) : 0);
               Duw_copy.coeffRef(id0,id0_connectB) = -cond.STEP/(2.*BC.valL);
            }
            else if (BC.typeL == string("Dirichlet")) Duw_copy.coeffRef(id0,id0) = 1.;
            if (BC.typeR == string("Neumann"))
            {
               Duw_copy.coeffRef(idN,idN) = 0.; // disconnect from the point itself
               Duw_copy.coeffRef(idN,idN_connectT) = cond.STEP/(2.*BC.valR);
               Duw_copy.coeffRef(idN,idN_connectB) = -cond.STEP/(2.*BC.valR);
            }
            else if (BC.typeR == string("Dirichlet")) Duw_copy.coeffRef(idN,idN) = 1.;

            // disconnect from default connections
            Duw_copy.coeffRef(id0,id0_disconnectB) = 0.;
            Duw_copy.coeffRef(id0,id0_disconnectT) = 0.;
            Duw_copy.coeffRef(idN,idN_disconnectB) = 0.;
            Duw_copy.coeffRef(idN,idN_disconnectT) = 0.;

            // If we know the VALUE at the point, add the index of it of the guess/solution
            //   vector that we already know, i.e. we don't need to update them
            if (!update && BC.typeL == string("Dirichlet")) no_update.push_back(ID(size,0,cond.SIZEu,n_v,op_elem_num));
            if (!update && BC.typeR == string("Dirichlet")) no_update.push_back(ID(size,cond.SIZEu-1,cond.SIZEu,n_v,op_elem_num));
         }

         for (int n_u = 1; n_u < cond.SIZEu-1; n_u++) // loop through the top and bottom boundary points of the mesh
         {
            // indexes for the bottom side
            int id0 =             ID(size, n_u,   cond.SIZEu, 0, 0), // the index of the mesh point we're at
                id0_connectL =    ID(size, n_u-1, cond.SIZEu, 0, 0), // index of the left point to connect to
                id0_connectR =    ID(size, n_u+1, cond.SIZEu, 0, 0), // index of the right point to connect to
                // we need to disconnect from the points that the default D matrix has
                id0_disconnectL = ID(size, n_u-1, cond.SIZEu, 1, 0), // index of the left point to disconnect from
                id0_disconnectR = ID(size, n_u+1, cond.SIZEu, 1, 0); // index of the right point to disconnect from
            
            // indexes for the top side
            int idN =             ID(size, n_u,   cond.SIZEu, cond.SIZEv-1, 0), // the index of the mesh point we're at
                idN_connectL =    ID(size, n_u-1, cond.SIZEu, cond.SIZEv-1, 0), // index of the left point to connect to
                idN_connectR =    ID(size, n_u+1, cond.SIZEu, cond.SIZEv-1, 0), // index of the right point to connect to
                // we need to disconnect from the points that the default D matrix has
                idN_disconnectL = ID(size, n_u-1, cond.SIZEu, cond.SIZEv-2, 0), // index of the left point to disconnect from
                idN_disconnectR = ID(size, n_u+1, cond.SIZEu, cond.SIZEv-2, 0); // index of the right point to disconnect from

            // set the values at these indexes using the ghost points
            if (BC.typeB == string("Neumann"))
            {
               Duw_copy.coeffRef(id0,id0) = 0.; // disconnect from the point itself
               Duw_copy.coeffRef(id0,id0_connectR) = cond.STEP/(2.*BC.valB);
               Duw_copy.coeffRef(id0,id0_connectL) = -cond.STEP/(2.*BC.valB);
            }
            else if (BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id0,id0) = 1.;
            if (BC.typeT == string("Neumann"))
            {
               Duw_copy.coeffRef(idN,idN) = 0.; // disconnect from the point itself
               Duw_copy.coeffRef(idN,idN_connectR) = cond.STEP/(2.*BC.valB);
               Duw_copy.coeffRef(idN,idN_connectL) = -cond.STEP/(2.*BC.valB);
            }
            else if (BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(idN,idN) = 1.;

            // disconnect from default connections
            Duw_copy.coeffRef(id0,id0_disconnectL) = 0.;
            Duw_copy.coeffRef(id0,id0_disconnectR) = 0.;
            Duw_copy.coeffRef(idN,idN_disconnectL) = 0.;
            Duw_copy.coeffRef(idN,idN_disconnectR) = 0.;

            // If we know the VALUE at the point, add the index of it of the guess/solution
            //   vector that we already know, i.e. we don't need to update them
            if (!update && BC.typeB == string("Dirichlet")) no_update.push_back(ID(size,n_u,cond.SIZEu,0,op_elem_num));
            if (!update && BC.typeT == string("Dirichlet")) no_update.push_back(ID(size,n_u,cond.SIZEu,cond.SIZEv-1,op_elem_num));
         }

         // special case for the corners
         int id, id_disconnect;

         // Top left
         id = ID(size,0,cond.SIZEu,cond.SIZEv-1,0);
         id_disconnect = ID(size,1,cond.SIZEu,cond.SIZEv-2,0);
         Duw_copy.coeffRef(id,id_disconnect) = 0.;
              if (BC.typeL == string("Neumann") && BC.typeT == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.valL*BC.valB);
         else if (BC.typeL == string("Dirichlet") || BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.;
               // TODO: determine if this assumption is correct, that
               //    if the function value is given for one side, we
               //    don't have to worry about the derivative condition.
               //    (x4, below)
         
         // Top right
         id = ID(size,cond.SIZEu-1,cond.SIZEu,cond.SIZEv-1,0);
         id_disconnect = ID(size,cond.SIZEu-2,cond.SIZEu,cond.SIZEv-2,0);
         Duw_copy.coeffRef(id,id_disconnect) = 0.;
              if (BC.typeR == string("Neumann") && BC.typeT == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.valR*BC.valB);
         else if (BC.typeR == string("Dirichlet") || BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; // here...

         // Bottom left
         id = ID(size,0,cond.SIZEu,0,0);
         id_disconnect = ID(size,1,cond.SIZEu,1,0);
         Duw_copy.coeffRef(id,id_disconnect) = 0.;
              if (BC.typeL == string("Neumann") && BC.typeB == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.valL*BC.valB);
         else if (BC.typeL == string("Dirichlet") || BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; //...here...

         // Bottom right
         id = ID(size,cond.SIZEu-1,cond.SIZEu,0,0);
         id_disconnect = ID(size,cond.SIZEu-2,cond.SIZEu,1,0);
         Duw_copy.coeffRef(id,id_disconnect) = 0.;
              if (BC.typeR == string("Neumann") && BC.typeB == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.valR*BC.valB);
         else if (BC.typeR == string("Dirichlet") || BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; //...and here
      
         // // If we know the VALUE at the point, add the index of it of the guess/solution
         // //   vector that we already know, i.e. we don't need to update them
         // if (BC.typeL == string("Dirichlet") || BC.typeT == string("Dirichlet")) no_update.push_back(ID(sz,0,        cond.SIZEu,cond.SIZEv-1,op_elem_num));
         // if (BC.typeR == string("Dirichlet") || BC.typeT == string("Dirichlet")) no_update.push_back(ID(sz,size[0]-1,cond.SIZEu,cond.SIZEv-1,op_elem_num));
         // if (BC.typeL == string("Dirichlet") || BC.typeB == string("Dirichlet")) no_update.push_back(ID(sz,0,        cond.SIZEu,0,        op_elem_num));
         // if (BC.typeR == string("Dirichlet") || BC.typeB == string("Dirichlet")) no_update.push_back(ID(sz,size[0]-1,cond.SIZEu,0,        op_elem_num));

         // cout << "\t\t\tDuw_BD: success" << endl;
         return Duw_copy;
      }
   
      // WRITE ADDITIONAL 'D**_BD' METHODS HERE

      void BuildProblem() {
         Build_D_Matrices();
         // cout << "build matrices: done" << endl;
         op_vector.resize(size*OP_size); // initialize OP vetor
         // cout << "op_vector resize: done" << endl;
         SolverMatrix = BuildSolverMatrix();
         // cout << "solver matrix: done" << endl;
      }

      // use the relaxation method and Anderson Acceleration to solve
      void Solve(Matrix<Scalar_type,-1,1>& guess) {
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
         else if (solver.info() == Eigen::NumericalIssue) { // for debugging non-invertable matrices
            cout << "Solver: numerical issues" << endl;
            // for (int k=0; k < SolverMatrix.outerSize(); ++k) for (SparseMatrix<Scalar_type>::InnerIterator it(SolverMatrix,k); it; ++it) cout << "(" << it.row() << "," << it.col() << ")\t";
            return;
         }

         // time to prepare solving method
         auto end = std::chrono::system_clock::now();
         auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
         cout << "\ttime: " << elapsed.count() << " seconds." << endl << endl;

         // update_OP_Vector(f); // prepare the matrix for RHS

         int cts = 0; // count loops
         double err;  // to store current error
         Matrix<Scalar_type,-1,1> rhs = RHS(); // the right hand side

         // the acceleration object
         converg_acceler<Matrix<Scalar_type,-1,1>> Con_Acc(no_update); // Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update)
         
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

            op_vector = f; // update the matrix for RHS
            rhs = RHS();   // update rhs
            cts++;         // increment counter

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
      //             __x__|__y__|_Eta_uu_|_Eta_uv_| ...
      //              ... | ... |  ...   |   ...  | ...
      // separates the components from solution...storing real and imag parts ==> up to 18
      void WriteToFile(Matrix<Scalar_type,-1,1>& vec, string file_name) {
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

      void WriteToFile_single_vector(VectorXd& vec, string file_name) {
         std::ofstream data (file_name);
         if (data.is_open()) {
            for (int i = 0; i < vec.size(); i++) {
               string line = to_string(i*cond.STEP) + string("\t") + to_string(vec(i));
               data << line << endl;
            }
         }
         else cout << "Unable to open file: " << file_name << endl;
      }

      void WriteToFile_single_vector(VectorXcd& vec, string file_name) {
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
// ==================================================================

// Now include the real and comples class header files,
//   since the parent class are defined
// #include "real_classes.hpp"
#include "complex_classes.hpp"

#endif