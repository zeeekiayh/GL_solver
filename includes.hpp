#ifndef _includes
#define _includes

#include <iostream>
#include <fstream>
#include <string>
#include "math.h"
#include <complex>
#include <vector>
#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>
#include "ConvergenceAccelerator.hpp"

using std::cout;
using std::endl;
using std::cin;
using std::complex;
using std::ostream;
using std::vector;
using std::string;
using namespace Eigen;
using namespace std::complex_literals; // for easy comlpex notation

typedef Triplet<complex<double>> Tr;

// constants used in the GL equations
struct GL_param { double K1, K2, K3, B1, B2, B3, B4, B5, A; };

// values used for set up, initial conditions, and algorithmic parameters
struct in_conditions
{
   int SIZEu,    // the number of points...
       SIZEv,    // "" ...
       SIZEw,    // "" in orthogonal directions
       N_loop,   // count limit for iterations
       maxStore, // max num vectors kept in accerleration
       wait;     // iterations to wait before using acceleration
   double ACCUR, // the minimum desired accuracy
          STEP,  // step size
          rel_p, // relaxation param
          bLeft, // slip length on left side
          bRight,//   "   "   on right side
          bTop,  //   "   "   on top side
          bBott; //   "   "   on bottom side
   // string BC0;   // boundary conditions at 0: "specular" or "diffuse"
};

void Write_To_File(VectorXcd& f, string f_name_real, string f_name_imag)
{
   // the real part
   std::ofstream data_r (f_name_real);
   if (data_r.is_open())
   {
      for (int i = 0; i < f.size(); i++) data_r << f(i).real() << "\n";
      data_r.close();
   }
   else cout << "Unable to open file: " << f_name_real << endl;

   // the imaginary part
   std::ofstream data_c (f_name_imag);
   if (data_c.is_open())
   {
      for (int i = 0; i < f.size(); i++) data_c << f(i).imag() << "\n";
      data_c.close();
   }
   else cout << "Unable to open file: " << f_name_imag << endl;
}

// returns the unique id corresponding to each op-component in the mesh.
//    n_u: mesh index 1
//    n_v: mesh index 2
//    i:   vector (OP) index
//    size: num elements in mesh
int ID(int size, int n_u, int n_u_max, int n_v, int i) { return size*i + n_u_max*n_v + n_u; }

// ==================================
   template <typename Container_type>
   struct OrderParam
   {
      // Container_type can be something like complex<double> or Vector3cd
      // If your order parameter is a matrix, it must be flattened as a Vector3cd

      // store it as a vector; the OP at one point
      int num_comp = 1; // we'll assume it's 1, unless in the derived class
      Container_type OP;

      OrderParam() {}
      OrderParam(int n): num_comp(n) {  }
      void Set_OP(Container_type op) { OP = op; }
      void initialize(int n)
      {
         num_comp = n;
         if (num_comp > 1) OP.resize(num_comp);
      }
      complex<double>& operator() (int i)
      {
         if (i > num_comp-1) // check the range first
         {
            cout << "ERROR: index out of range OrderParam::operator()" << endl;
            // return 0.;
         }
         // will these all have to be specialized??
         // if (num_comp == 1) return OP; // if it's just one-component OP, return that value
         else return OP(i);            // otherwise, it is a vector, so return it's component
      }
   };

   template<class Container_type>
   ostream& operator<< (ostream& os, const OrderParam<Container_type>& OP)
   {
      for (int i = 0; i < OP.num_comp; i++)
      {
         os << OP.OP(i);
         if (i+1 < OP.num_comp) os << "\t";
      }
      return os;
   }
// ==================================


// ======================================================
// to allow for the specialized struct, define a template
   template<typename Container_type>
   struct MultiComponentOrderParam : public OrderParam<Container_type> {};

// a structure specialized for multi-component order parameters
   template<>
   struct MultiComponentOrderParam<VectorXcd> : public OrderParam<VectorXcd>
   {
      MultiComponentOrderParam(int n) { initialize(n); }
      void Set_OP(Matrix3cd op) // will we actually use this??
      {
         // flatten the matrix (row major)
         int i = 0; // count the values put into the vector OP

         for (int a = 0; a < 3; a++) {
            for (int j = 0; j < 3; j++) {
               if (abs(op(a,j)) > pow(10,-8)) { // if not 0
                  if (i > num_comp) cout << "WARNING: more elements in matrix than specified by num_comp." << endl;
                  else {
                     this->OP(i) = op(a,j);
                     i++;
                  }
               }
            }
         } // for's
      }
      Matrix3cd GetMatrixForm_He3Defect() // this function is specific to one OP structure
      {
         Matrix3cd A(3,3);
         A << OP(0), 0.,    OP(1),
              0.,    OP(2), 0.,
              OP(3), 0.,    OP(4);
         return A;
      }
   };
// ======================================================


// ===========================================
// Boundary condition structure for a sigle OP
   struct Bound_Cond
   {
      string typeB, typeT, typeL, typeR; // type of BC: Dirichlet (value) or Neumann (derivative)
      double bB, bT, bL, bR; // for Neumann BC:   slip length
                             // for Dirichlet BC: function value
   };
// ===========================================


// ================================================================
// a class to hold the mesh we use in FEM, which holds all the OP's
   template <typename Container_type>
   class OP_Matrix
   {
      public:
      OP_Matrix() {}
      // for multi-component OP's
      OP_Matrix(int n, int nRows, int nCols) { initialize(n, nRows, nCols); }
      // OP_Matrix(int n, int nRows, int nCols, int s3) {}

      void initialize(int n, int nRows, int nCols)
      {
         cout << "in OP_Martix::initialize()" << endl;
         size[0] = nRows; // initialize size variables
         size[1] = nCols;
         OP_size = n; // number of OP components
         cout << "member variables initialized" << endl;

         matrix.resize(size[0],size[1]); // initialize matrix
         vector.resize(size[0]*size[1]*OP_size); // initialize vector, for the whole thing (size = num_of_mesh_points * num_OP_components)
         cout << "this->matrix and this->vector resized" << endl;

         // initialize elements in 'matrix'
         for (int i = 0; i < size[0]; i++) {
            for (int j = 0; j < size[1]; j++) {
               matrix(i,j).initialize(OP_size);
            }
         }
         cout << "all matrix elements successfully initialized" << endl;

         // make the vector form available
         setVectorForm();
         cout << "made vector form available" << endl;
      }

      // derivative matrix methods
      SparseMatrix<complex<double>> Du2_BD(SparseMatrix<complex<double>>& Du2, double h, Bound_Cond BC)
      {
         // the matrix that we will edit and return to not modify the original
         SparseMatrix<complex<double>> Du2_copy = Du2;
         
         // loop through just the left and right boundary points of the mesh
         for (int n_v = 0; n_v < size[1]; n_v++)
         {
            // indexes for the left side
            int id0 =         ID(size[0]*size[1], 0, size[0], n_v, 0),
                id0_connect = ID(size[0]*size[1], 1, size[0], n_v, 0);
            
            // indexes for the right side
            int idN =         ID(size[0]*size[1], size[0]-1, size[0], n_v, 0),
                idN_connect = ID(size[0]*size[1], size[0]-2, size[0], n_v, 0);

            // set the values at these indexes using the ghost points
            if (BC.typeL == string("Neumann"))
            {
               Du2_copy.coeffRef(id0,id0) = -2. -2.*h/BC.bL;
               Du2_copy.coeffRef(id0,id0_connect) = 2.;
            }
            else if (BC.typeL == string("Dirichlet"))
            {
               cout << "";
               Du2_copy.coeffRef(id0,id0) = 1.;
            }

            if (BC.typeR == string("Neumann"))
            {
               Du2_copy.coeffRef(idN,idN) = -2. -2.*h/BC.bR; // TODO: is it -2h... or +2h... ?
               Du2_copy.coeffRef(idN,idN_connect) = 2.;
            }
            else if (BC.typeR == string("Dirichlet"))
            {
               cout << "";
               Du2_copy.coeffRef(idN,idN) = 1.;
            }
         }

         return Du2_copy;
      }
      
      SparseMatrix<complex<double>> Dv2_BD(SparseMatrix<complex<double>>& Dv2, double h, Bound_Cond BC)
      {
         // the matrix that we will edit and return to not modify the original
         SparseMatrix<complex<double>> Dv2_copy = Dv2;

         // loop through just the top and bottom boundary points of the mesh
         for (int n_u = 0; n_u < size[0]; n_u++)
         {
            // indexes for the bottom side
            int id0 =         ID(size[0]*size[1], n_u, size[0], 0, 0),
                id0_connect = ID(size[0]*size[1], n_u, size[0], 1, 0);
            
            // indexes for the top side
            int idN =         ID(size[0]*size[1], n_u, size[0], size[1]-1, 0),
                idN_connect = ID(size[0]*size[1], n_u, size[0], size[1]-2, 0);

            // set the values at these indexes using the ghost points
            if (BC.typeB == string("Neumann"))
            {
               Dv2_copy.coeffRef(id0,id0) = -2. -2.*h/BC.bB;
               Dv2_copy.coeffRef(id0,id0_connect) = 2.;
            }
            else if (BC.typeB == string("Dirichlet")) Dv2_copy.coeffRef(id0,id0) = 1.;

            if (BC.typeT == string("Neumann"))
            {
               Dv2_copy.coeffRef(idN,idN) = -2. -2.*h/BC.bT; // TODO: is it -2h... or +2h... ?
               Dv2_copy.coeffRef(idN,idN_connect) = 2.;
            }
            else if (BC.typeT == string("Dirichlet")) Dv2_copy.coeffRef(idN,idN) = 1.;
         }

         return Dv2_copy;
      }
      
      SparseMatrix<complex<double>> Duv_BD(SparseMatrix<complex<double>>& Duv, double h, Bound_Cond BC)
      {
         cout << "in Duv_BD()" << endl;
         // the matrix that we will edit and return to not modify the original
         SparseMatrix<complex<double>> Duv_copy = Duv;

         int sz = size[0]*size[1]; // size of D matrix (num of mesh points)

         // loop through just the boundary points of the mesh
         for (int n_v = 1; n_v < size[1]-1; n_v++) // loop through the left and right boundary points of the mesh
         {
            // indexes for the left side
            int id0 =             ID(sz, 0, size[0], n_v,   0), // the index of the mesh point we're at
                id0_connectB =    ID(sz, 0, size[0], n_v-1, 0), // index of the bottom point to connect to
                id0_connectT =    ID(sz, 0, size[0], n_v+1, 0), // index of the top point to connect to
               // we need to disconnect from the points that the default D matrix has
                id0_disconnectB = ID(sz, 1, size[0], n_v-1, 0), // index of the bottom point to disconnect from
                id0_disconnectT = ID(sz, 1, size[0], n_v+1, 0); // index of the top point to disconnect from
            
            // indexes for the right side
            int idN =             ID(sz, size[0]-1, size[0], n_v,   0), // the index of the mesh point we're at
                idN_connectB =    ID(sz, size[0]-1, size[0], n_v-1, 0), // index of the bottom point to connect to
                idN_connectT =    ID(sz, size[0]-1, size[0], n_v+1, 0), // index of the top point to connect to
               // we need to disconnect from the points that the default D matrix has
                idN_disconnectB = ID(sz, size[0]-2, size[0], n_v-1, 0), // index of the bottom point to disconnect from
                idN_disconnectT = ID(sz, size[0]-2, size[0], n_v+1, 0); // index of the top point to disconnect from

            // set the values at these indexes using the ghost points
            if (BC.typeL == string("Neumann"))
            {
               Duv_copy.coeffRef(id0,id0) = 0.; // disconnect from the point itself
               Duv_copy.coeffRef(id0,id0_connectT) = h/(2.*BC.bL);
               Duv_copy.coeffRef(id0,id0_connectB) = -h/(2.*BC.bL);
            }
            else if (BC.typeL == string("Dirichlet")) Duv_copy.coeffRef(id0,id0) = 1.;

            if (BC.typeR == string("Neumann"))
            {
               Duv_copy.coeffRef(idN,idN) = 0.; // disconnect from the point itself
               Duv_copy.coeffRef(idN,idN_connectT) = h/(2.*BC.bR);
               Duv_copy.coeffRef(idN,idN_connectB) = -h/(2.*BC.bR);
            }
            else if (BC.typeR == string("Dirichlet")) Duv_copy.coeffRef(idN,idN) = 1.;

            // disconnect from default connections
            Duv_copy.coeffRef(id0,id0_disconnectB) = 0.;
            Duv_copy.coeffRef(id0,id0_disconnectT) = 0.;
            Duv_copy.coeffRef(idN,idN_disconnectB) = 0.;
            Duv_copy.coeffRef(idN,idN_disconnectT) = 0.;
         }
         for (int n_u = 1; n_u < size[0]-1; n_u++) // loop through the top and bottom boundary points of the mesh
         {
            // indexes for the bottom side
            int id0 =             ID(sz, n_u,   size[0], 0, 0), // the index of the mesh point we're at
                id0_connectL =    ID(sz, n_u-1, size[0], 0, 0), // index of the left point to connect to
                id0_connectR =    ID(sz, n_u+1, size[0], 0, 0), // index of the right point to connect to
               // we need to disconnect from the points that the default D matrix has
                id0_disconnectL = ID(sz, n_u-1, size[0], 1, 0), // index of the left point to disconnect from
                id0_disconnectR = ID(sz, n_u+1, size[0], 1, 0); // index of the right point to disconnect from
            
            // indexes for the top side
            int idN =             ID(sz, n_u,   size[0], size[1]-1, 0), // the index of the mesh point we're at
                idN_connectL =    ID(sz, n_u-1, size[0], size[1]-1, 0), // index of the left point to connect to
                idN_connectR =    ID(sz, n_u+1, size[0], size[1]-1, 0), // index of the right point to connect to
               // we need to disconnect from the points that the default D matrix has
                idN_disconnectL = ID(sz, n_u-1, size[0], size[1]-2, 0), // index of the left point to disconnect from
                idN_disconnectR = ID(sz, n_u+1, size[0], size[1]-2, 0); // index of the right point to disconnect from

            // set the values at these indexes using the ghost points
            if (BC.typeB == string("Neumann"))
            {
               Duv_copy.coeffRef(id0,id0) = 0.; // disconnect from the point itself
               Duv_copy.coeffRef(id0,id0_connectR) = h/(2.*BC.bB);
               Duv_copy.coeffRef(id0,id0_connectL) = -h/(2.*BC.bB);
            }
            else if (BC.typeB == string("Dirichlet")) Duv_copy.coeffRef(id0,id0) = 1.;

            if (BC.typeT == string("Neumann"))
            {
               Duv_copy.coeffRef(idN,idN) = 0.; // disconnect from the point itself
               Duv_copy.coeffRef(idN,idN_connectR) = h/(2.*BC.bT);
               Duv_copy.coeffRef(idN,idN_connectL) = -h/(2.*BC.bT);
            }
            else if (BC.typeT == string("Dirichlet")) Duv_copy.coeffRef(idN,idN) = 1.;

            // disconnect from default connections
            Duv_copy.coeffRef(id0,id0_disconnectL) = 0.;
            Duv_copy.coeffRef(id0,id0_disconnectR) = 0.;
            Duv_copy.coeffRef(idN,idN_disconnectL) = 0.;
            Duv_copy.coeffRef(idN,idN_disconnectR) = 0.;
         }

         cout << "finished for-loops" << endl;

         // special case for the corners
         int id, id_disconnect;

         // Top left
         id_disconnect = ID(sz,1,size[0],size[1]-2,0);
         id = ID(sz,0,size[0],size[1]-1,0);
         Duv_copy.coeffRef(id,id_disconnect) = 0.;
              if (BC.typeL == string("Neumann") && BC.typeT == string("Neumann"))     Duv_copy.coeffRef(id,id) = h*h/(BC.bL*BC.bT);
         else if (BC.typeL == string("Dirichlet") || BC.typeT == string("Dirichlet")) Duv_copy.coeffRef(id,id) = 1.;
               // TODO: determine if this assumption is correct, that
               //    if the function value is given for one side, we
               //    don't have to worry about the derivative condition.
               //    (x4, below)
         
         cout << "here 1" << endl;

         // Top right
         id_disconnect = ID(sz,size[0]-2,size[0],size[1]-2,0);
         id = ID(sz,size[0]-1,size[0],size[1]-1,0);
         Duv_copy.coeffRef(id,id_disconnect) = 0.;
              if (BC.typeR == string("Neumann") && BC.typeT == string("Neumann"))     Duv_copy.coeffRef(id,id) = h*h/(BC.bR*BC.bT);
         else if (BC.typeR == string("Dirichlet") || BC.typeT == string("Dirichlet")) Duv_copy.coeffRef(id,id) = 1.; // here...

         cout << "here 2" << endl;

         // Bottom left
         id_disconnect = ID(sz,1,size[0],1,0);
         id = ID(sz,0,size[0],0,0);
         Duv_copy.coeffRef(id,id_disconnect) = 0.;
              if (BC.typeL == string("Neumann") && BC.typeB == string("Neumann"))     Duv_copy.coeffRef(id,id) = h*h/(BC.bL*BC.bB);
         else if (BC.typeL == string("Dirichlet") || BC.typeB == string("Dirichlet")) Duv_copy.coeffRef(id,id) = 1.; //...here...

         cout << "here 3" << endl;

         // Bottom right
         id_disconnect = ID(sz,size[0]-2,size[0],1,0);
         id = ID(sz,size[0]-1,size[0],0,0);
         Duv_copy.coeffRef(id,id_disconnect) = 0.;
              if (BC.typeR == string("Neumann") && BC.typeB == string("Neumann"))     Duv_copy.coeffRef(id,id) = h*h/(BC.bR*BC.bB);
         else if (BC.typeR == string("Dirichlet") || BC.typeB == string("Dirichlet")) Duv_copy.coeffRef(id,id) = 1.; //...and here

         cout << "here 4" << endl;

         return Duv_copy;
      }
      
      // Convert the OP matrix into a vector
      // This functions may need to be used repeatedly
      //   as the guess is updated...or we could just
      //   use the &getVector() method to change it directly
      void setVectorForm()
      {
         for (int vi = 0; vi < OP_size; vi++) {
            for (int row = 0; row < size[0]; row++) {
               for (int col = 0; col < size[1]; col++) {
                  vector(ID(size[0]*size[1],row,size[0],col,vi)) = matrix(row,col)(vi);
               }
            }
         }
      }

      VectorXcd getVector() const { return vector; }
      VectorXcd& getVector()      { return vector; }

      // free-energy functions
      complex<double> f_bulk();
      complex<double> f_grad();
      double free_energy();

      protected:
      int OP_size; // number of components in a single OP
      int size[2]; // to hold the 2 sizes, in each orthogonal direction along the mesh
      // int size[3]; // to hold the possible 3 sizes

      // to hold the OP at each point on the mesh
      Matrix<OrderParam<Container_type>,-1,-1> matrix;

      // the vector form of matrix
      VectorXcd vector;
   };

   template <typename Container_type>
   class MultiComponent_OP_Matrix: public OP_Matrix<Container_type> {};
   template <>
   class MultiComponent_OP_Matrix<VectorXcd>: public OP_Matrix<VectorXcd>
   {
      public:
      // User-defined methods to build the solver matrix and the rhs vector
      SparseMatrix<complex<double>> SolverMatrix_He3Defect(GL_param gl,
                                                           SparseMatrix<complex<double>>& Du2,
                                                           SparseMatrix<complex<double>>& Dv2,
                                                           SparseMatrix<complex<double>>& Duv,
                                                           double h, Bound_Cond Axx,
                                                           Bound_Cond Axz, Bound_Cond Ayy,
                                                           Bound_Cond Azx, Bound_Cond Azz)
      {
         cout << "building SolverMatrix_He3Defect" << endl;
         // to make the code cleaner, define some constants
         double K123 = gl.K1+gl.K2+gl.K3,
                K23  = gl.K2+gl.K3;

         // define matrices to use, with their sizes
         MatrixXcd SolverMatrix(Du2.rows()*5,Du2.cols()*5); // the matrix to be used by the solver
         MatrixXcd zero(Du2.rows(),Du2.cols()); // a zero matrix to fill in the spaces with the comma initializer
         zero.setZero(); // make sure it's zero
         cout << "matrices prepared" << endl;

         // define each non-zero 'element'
         MatrixXcd elem_00 = MatrixXcd(K123*Du2_BD(Du2,h,Axx)+gl.K1*Dv2_BD(Dv2,h,Axx));
            cout << "elem_00 initialized" << endl;
         MatrixXcd elem_10 = MatrixXcd(K23*Duv_BD(Duv,h,Axx));
            cout << "elem_10 initialized" << endl;
         MatrixXcd elem_01 = MatrixXcd(K23*Duv_BD(Duv,h,Axz));
            cout << "elem_01 initialized" << endl;
         MatrixXcd elem_11 = MatrixXcd(K123*Duv_BD(Duv,h,Axz)+gl.K1*Du2_BD(Du2,h,Axz));
            cout << "elem_11 initialized" << endl;
         MatrixXcd elem_22 = MatrixXcd(gl.K1*(Du2_BD(Du2,h,Ayy)+Dv2_BD(Dv2,h,Ayy)));
            cout << "elem_22 initialized" << endl;
         MatrixXcd elem_33 = MatrixXcd(K123*Du2_BD(Du2,h,Azx)+gl.K1*Dv2_BD(Dv2,h,Azx));
            cout << "elem_33 initialized" << endl;
         MatrixXcd elem_43 = MatrixXcd(K23*Duv_BD(Duv,h,Azx));
            cout << "elem_43 initialized" << endl;
         MatrixXcd elem_34 = MatrixXcd(K23*Duv_BD(Duv,h,Azz));
            cout << "elem_34 initialized" << endl;
         MatrixXcd elem_44 = MatrixXcd(K123*Duv_BD(Duv,h,Azz)+gl.K1*Du2_BD(Du2,h,Azz));
         cout << "matrix elements defined" << endl;

         // use the comma initializer to build the matrix
         SolverMatrix << elem_00, elem_10, zero,    zero,    zero,
                         elem_10, elem_11, zero,    zero,    zero,
                         zero,    zero,    elem_22, zero,    zero,
                         zero,    zero,    zero,    elem_33, elem_34,                          
                         zero,    zero,    zero,    elem_43, elem_44;
         cout << "comma initializer success for SolverMatrix" << endl;

         // turn the matrix into a sparse matrix
         return SolverMatrix.sparseView(1,pow(10,-8));
      }

      VectorXcd RHS_He3Defect(GL_param gl)
      {
         Matrix3cd A, A_T, A_dag, A_conj;
         VectorXcd rhs(vector.size());

         // loop through all the OP components in the mesh
         for (int vi = 0; vi < OP_size; vi++) {
            for (int row = 0; row < size[0]; row++) {
               for (int col = 0; col < size[1]; col++) {
                  complex<double> A_mui = vector(ID(size[0]*size[1],row,size[0],col,vi));
                  
                  A = matrix(row,col).GetMatrixForm_He3Defect();
                  A_T = A.transpose();
                  A_dag = A.adjoint();
                  A_conj = A.conjugate();

                  complex<double> val = 2*gl.B1*(A*A_T).trace()*conj(A_mui)
                                       +2*gl.B2*(A*A_dag).trace()*A_mui
                                       +2*gl.B3*(A*A_T*A_conj)(floor(2*vi/3),(2*vi)%3)
                                       +2*gl.B4*(A*A_dag*A)(floor(2*vi/3),(2*vi)%3)
                                       +2*gl.B5*(A_conj*A_T*A)(floor(2*vi/3),(2*vi)%3)
                                       +gl.A*A_mui;
                  // Where (floor(2*vi/3),(2*vi)%3) gives us (row#, col#) in the
                  //   OP matrix, given a value for vi: the # of element in the
                  //   vector. This is specific to this one case with 5 OP components
                  //   that are arranged in the corners and center.

                  // add val to rhs, in the matching location
                  rhs(ID(size[0]*size[1],row,size[0],col,vi)) = val;
               }
            }
         }
         return rhs;
      }

      private:
      Matrix<MultiComponentOrderParam<VectorXcd>,-1,-1> matrix;
   };
// ================================================================


// =============================================
// specialized classes for single-component OP's
   template <>
   class OP_Matrix<complex<double>>
   {
      //
   };
   template <>
   class OP_Matrix<double>
   {
      //
   };
// =============================================

// TODO: modify the grad terms to calculate it based on the mesh
// Calculate the free-energy loss by the integral of the energy density
// This value has been normalized because the deltas were calculated as
//    normalized. Dividing f by (alpha(T) * delta_0^2) and simplifying.
// Gives a warning if the result has ~non-zero imaginary part.
// double Free_energy()
// {
//    if (!solution.size())
//    {
//       cout << "ERROR: cannot calculate free-energy without a solution." << endl;
//       return 0.;
//    }
//    complex<double> I = 0; // start the integral sum at 0
//    complex<double> f_bulk, f_bulk_prev = 0.;
//    complex<double> f_grad, f_grad_prev = 0.;
//    VectorXcd integ(size-2); // value of the integral over distance--to plot
//    // calculate the first step
//    f_bulk = F_Bulk(0);
//    f_grad = F_Grad(0,1);
//    I += ( f_bulk + f_grad - 1. )*cond.STEP;
//    for (int i = 1; i <= size-2; i++)
//    {
//       // set the previous values
//       f_bulk_prev = f_bulk;
//       f_grad_prev = f_grad;
//       f_bulk = F_Bulk(i);
//       // calculate the gradient term
//       f_grad = F_Grad(i-1,i+1);
//       // use a rectangular integral approximation, centered at the midpoints
//       I += ( (f_bulk+f_bulk_prev + f_grad+f_grad_prev)/2. - 1. )*cond.STEP;
//       integ(i-1) = (f_bulk+f_bulk_prev + f_grad+f_grad_prev)/2.;//I;
//    }
//    // calculate the last step
//    f_bulk = F_Bulk(size-1);
//    f_grad = F_Grad(size-2,size-1);
//    I += ( f_bulk + f_grad - 1. )*cond.STEP;
//    // save the integrand vector to plot and inspect
//    Write_To_File(integ,"integ_c.txt","integ_r.txt"); // using the non-member function
//    cout << "The final value of f/f0 = " << integ(integ.size()-1) << endl;
//    if (I.imag() >= pow(10,-8)) cout << "WARNING: imaginary part of the free-energy is not zero." << endl;
//    return I.real();
// }
// calculate the normalized bulk free-energy density
// complex<double> F_Bulk(int i)
// {
//    // calculate the used forms of A
//    Matrix<complex<double>,3,3> A = M_index(OP,i),
//                               AT = A.transpose(),
//                               A_dag = A.adjoint(),
//                               A_conj = A.conjugate();
//    double Beta_B = gl.B1+gl.B2 + (gl.B3+gl.B4+gl.B5)/3.;
//    return -( gl.B1*pow( abs((A * AT).trace()), 2)
//             +gl.B2*pow( (A * A_dag).trace(), 2)
//             +gl.B3*(A * AT * A_conj * A_dag).trace()
//             +gl.B4*(A * A_dag * A * A_dag).trace()
//             +gl.B5*(A * A_dag * A_conj * AT).trace()
//          )/(Beta_B*9.)
//          +2./3.*(A * A_dag).trace();
// }
// TODO: make it calculate the gradient using the central difference
//       derivative for each internal point on the mesh...do we pass
//       in the mesh? or just pass in the relating indexes?
// calculate the normalized gradient free-energy density
// complex<double> F_Grad(int m, int n)
// {
//    if (m < 0) { cout << "ERROR: F_Grad 'm' must be >= 0." << endl; return 0.0; }
//    complex<double> k1 = 0., k2 = 0., k3 = 0.; // grad term sums
//    // used center difference derivatives
//    Matrix<complex<double>,3,3> A_next = M_index(OP,n), A_prev = M_index(OP,m);
//    for (int a = 0; a < 3; a++)
//    {
//       for (int k = 0; k < 3; k++)
//       {
//          for (int j = 0; j < 3; j++)
//          {
//             // these derivatives are divided by their step size when used in Free_Energy()
//             if (var_mat(a,j)[k])                 k1 += ( conj(A_next(a,j) - A_prev(a,j)) ) * ( A_next(a,j) - A_prev(a,j) );
//             if (var_mat(a,k)[j]*var_mat(a,j)[k]) k2 += ( conj(A_next(a,k) - A_prev(a,k)) ) * ( A_next(a,j) - A_prev(a,j) );
//             if (var_mat(a,j)[k]*var_mat(a,k)[j]) k3 += ( conj(A_next(a,j) - A_prev(a,j)) ) * ( A_next(a,k) - A_prev(a,k) );
//          }
//       }
//    }
//    return -2./3.*(k1+k2+k3)/(pow(2*cond.STEP,2)); // divide by (2h)^2 becasue there
//                   //  is a product of 2 derivatives, but we have double step size
// }

// ===============================
   template <class Container_type>
   class GL_Solver
   {
      // * * * * * * *
      // * VARIABLES *
      // * * * * * * *
      private:
      OP_Matrix<Container_type> op_matrix; // the matrix of OP at each mesh point

      protected:
      int size; // the size of the D matrices and the OP-component vectors
      GL_param gl; // temperature-dependent parameters for the GL equation
      in_conditions cond; // struct of all the BC's and other parameters for the methods
      VectorXcd solution; // to store the solution to the GL equ. (in the single vector form)

      SparseMatrix<complex<double>> A; // solver matrix
      SparseMatrix<complex<double>> Du2, Dv2, Duv; // derivative matrices

      public:
      // * * * * * * * * * * * * * * * *
      // * CONSTRUCTORS & DECSTRUCTOR  *
      // * * * * * * * * * * * * * * * *
      GL_Solver() {};
      GL_Solver(string conditions_file)
      {
         ReadConditions(conditions_file);
         size = cond.SIZEu*cond.SIZEv;
      }
      // ~GL_Solver() {};

      // * * * * * *
      // * METHODS *
      // * * * * * *
      void Solve(VectorXcd&);
      void BuildProblem(int,Bound_Cond,Bound_Cond,Bound_Cond,Bound_Cond,Bound_Cond);

      // Write all components of the OP, all into one file, of the form:
      //             __x__|__y__|_Axx_|_Axy_| ...
      //              ... | ... | ... | ... | ...
      //   separates the components from solution...storing real and imag parts ==> up to 18 files
      void WriteToFile(string f_name_real, string f_name_imag)
      {
         // the real part
         std::ofstream data_r (f_name_real);
         if (data_r.is_open())
         {
            for (int i = 0; i < solution.size(); i++) data_r << solution(i).real() << "\n";
            data_r.close();
         }
         else cout << "Unable to open file: " << f_name_real << endl;

         // the imaginary part
         std::ofstream data_c (f_name_imag);
         if (data_c.is_open())
         {
            for (int i = 0; i < solution.size(); i++) data_c << solution(i).imag() << "\n";
            data_c.close();
         }
         else cout << "Unable to open file: " << f_name_imag << endl;
      }

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
                  else if (ls == string("ACCUR"))     conditions >> cond.ACCUR;
                  else if (ls == string("h"))         conditions >> cond.STEP;
                  else if (ls == string("num loops")) conditions >> cond.N_loop;
                  else if (ls == string("rel_p"))     conditions >> cond.rel_p;
                  else if (ls == string("maxStore"))  conditions >> cond.maxStore;
                  else if (ls == string("wait"))      conditions >> cond.wait;
                  else if (ls == string("bLeft"))     conditions >> cond.bLeft;
                  else if (ls == string("bRight"))    conditions >> cond.bRight;
                  else if (ls == string("bTop"))      conditions >> cond.bTop;
                  else if (ls == string("bBott"))     conditions >> cond.bBott;
                  else if (ls == string("betas"))
                  {
                     conditions >> gl.B1;
                     conditions >> gl.B2;
                     conditions >> gl.B3;
                     conditions >> gl.B4;
                     conditions >> gl.B5;
                  }
                  else if (ls == string("alpha")) conditions >> gl.A;
                  else if (ls == string("Ks"))
                  {
                     conditions >> gl.K1;
                     conditions >> gl.K2;
                     conditions >> gl.K3;
                  }
               }
            }
            conditions.close();
         }
         else cout << "Unable to open file:" << conditions_file << endl;
      }

      // Build the derivative matrices
      void Build_D_Matrices()
      {
         // make vectors to hold the triplets of coefficients
         vector<Tr> coeffs_x2, coeffs_z2, coeffs_xz; // coeffs_y2, coeffs_xy, coeffs_yz

         // u, v, w are orthogonal basis for the system
         //   the n's are indexes in each direction

         // Loop through the entire mesh
         for (int n_u = 0; n_u < cond.SIZEu; n_u++)
         {
            for (int n_v = 0; n_v < cond.SIZEv; n_v++)
            {
               int id = ID(size,n_u,cond.SIZEu,n_v,0);
               
               InsertCoeff_Du2(id, n_u,   n_v, -2., coeffs_x2);
               InsertCoeff_Du2(id, n_u-1, n_v,  1., coeffs_x2);
               InsertCoeff_Du2(id, n_u+1, n_v,  1., coeffs_x2);
               
               InsertCoeff_Dv2(id, n_u, n_v,  -2., coeffs_z2);
               InsertCoeff_Dv2(id, n_u, n_v-1, 1., coeffs_z2);
               InsertCoeff_Dv2(id, n_u, n_v+1, 1., coeffs_z2);
               
               InsertCoeff_Duv(id, n_u-1, n_v-1, 1./4., coeffs_xz);
               InsertCoeff_Duv(id, n_u+1, n_v-1,-1./4., coeffs_xz);
               InsertCoeff_Duv(id, n_u-1, n_v+1,-1./4., coeffs_xz);
               InsertCoeff_Duv(id, n_u+1, n_v+1, 1./4., coeffs_xz);
               
            }
         } // end for's

         // initialize the D's by size
         Du2.resize(size,size);
         Dv2.resize(size,size);
         Duv.resize(size,size);

         // build all the D's from their coefficient triplet-vectors
         Du2.setFromTriplets(coeffs_x2.begin(), coeffs_x2.end());
         Dv2.setFromTriplets(coeffs_z2.begin(), coeffs_z2.end());
         Duv.setFromTriplets(coeffs_xz.begin(), coeffs_xz.end());
      }

      // Insert method for the Dx^2 matrix derivatives
      void InsertCoeff_Du2(int id, int u, int v, complex<double> weight, vector<Tr>& coeffs)
      {
         // id is the index of the connecting element, so
         //   id1 is the index of the element being inserted
         int id1 = ID(size,u,cond.SIZEu,v,0);

         if (u == -1 || u == cond.SIZEu);  // would add boundary conditions here, but we'll use ghost points, so do nothing
         // we'll just insert default values, filling the spaces, but they will be modified later for each OP-component
         else coeffs.push_back(Tr(id,id1,weight));
      }

      // insert method for the Dz^2 matrix derivatives
      void InsertCoeff_Dv2(int id, int u, int v, complex<double> weight, vector<Tr>& coeffs)
      {
         // id is the index of the connecting element, so
         //   id1 is the index of the element being inserted
         int id1 = ID(size,u,cond.SIZEu,v,0);

         if (v == -1 || v == cond.SIZEv);  // would add boundary conditions here, but we'll use ghost points, so do nothing
         else coeffs.push_back(Tr(id,id1,weight));
      }

      // insert method for the mixed derivative matrices
      void InsertCoeff_Duv(int id, int u, int v, complex<double> weight, vector<Tr>& coeffs)
      {
         // id is the index of the connecting element, so
         //   id1 is the index of the element being inserted
         int id1 = ID(size,u,cond.SIZEu,v,0);

              if (u == -1 || u == cond.SIZEu); // would add boundary conditions here,
         else if (v == -1 || v == cond.SIZEv); //  but we'll use ghost points, so do nothing
         else coeffs.push_back(Tr(id,id1,weight));
      }

   }; // GL_solver class

// derived, multi-component GL solver class
   template<class Container_type>
   class MultiComponent_GL_Solver : public GL_Solver<Container_type> {};
   template<>
   class MultiComponent_GL_Solver<VectorXcd> : public GL_Solver<VectorXcd>
   {
      public:
      MultiComponent_GL_Solver(string conditions_file)
      {
         ReadConditions(conditions_file);
         size = cond.SIZEu*cond.SIZEv;
      }

      void BuildProblem(int n, Bound_Cond Axx, Bound_Cond Axz,
                        Bound_Cond Ayy, Bound_Cond Azx,
                        Bound_Cond Azz)
      {
         Build_D_Matrices();
         // cout << "Du2 =\n" << Du2.real() << endl;
         cout << "Dv2 =\n" << Dv2.real() << endl;
         // cout << "Duv =\n" << Duv.real() << endl;

         this->op_matrix.initialize(n,cond.SIZEu,cond.SIZEv);
         cout << "initialized op_matrix" << endl;
         A = op_matrix.SolverMatrix_He3Defect(gl,Du2,Dv2,Duv,cond.STEP,Axx,Axz,Ayy,Azx,Azz);
         cout << "got SolverMatrix_He3Defect" << endl;
      }

      // TODO...?
      // use the relaxation method and Anderson Acceleration to solve
      void Solve(VectorXcd& guess)
      {
         VectorXcd f = guess, df(guess.size()); // initialize vectors

         // TODO: change this...it will be different with all the
         //       op-components in the vector
         // the elements that we don't want changed in the acceleration method
         vector<int> no_update;
         no_update.push_back(1);
         no_update.push_back(guess.size()-2);
         no_update.push_back(guess.size()-1);

         // use LU decomposition to solve the system
         SparseLU<SparseMatrix<complex<double>>, COLAMDOrdering<int> > solver;
         solver.analyzePattern(A); // TODO: check to see if needed
         solver.factorize(A);

         // check to see if the solver failed
         if (solver.info() == Eigen::Success) cout << "Solver: success" << endl;
         else if (solver.info() == Eigen::NumericalIssue) cout << "Solver: numerical issues" << endl;

         // loop until f converges or until it's gone too long
         int cts = 0; // count loops
         double err;  // to store current error
         VectorXcd rhs = op_matrix.RHS_He3Defect(gl); // the right hand side

         // the acceleration object
         converg_acceler<VectorXcd> Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update);

         do { // use relaxation
            df = solver.solve(rhs)-f; // find the change in f

            // use Anderson Acceleration to converge faster
            Con_Acc.next_vector<Matrix<dcomplex,-1,-1>>(f,df,err);

            rhs = op_matrix.RHS_He3Defect(gl); // update rhs
            cts++;
         } while(err > cond.ACCUR && cts < cond.N_loop);

         if (err < cond.ACCUR) cout << "Found solution:" << endl;
         cout << "   iterations = " << cts << endl;
         cout << "   relative error = " << err << endl;

         solution = f;
         // return f; // return the last value for f
      }

      private:
      MultiComponent_OP_Matrix<VectorXcd> op_matrix;
   };
// ===============================

#endif
