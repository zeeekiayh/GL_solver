#ifndef structures_hpp_
#define structures_hpp_

#include <iostream>
#include <iomanip> // for std::setprecision
#include <fstream> // for file in/out
#include <math.h>
#include <complex> 
#include <typeinfo>

#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>

#include "Kstorage.hpp"
// ===========================================
// type definitions for complex OP
   typedef std::complex<double> T_scalar;
   typedef Eigen::VectorXcd  T_vector;
   typedef Eigen::MatrixXcd  T_matrix;
   typedef Eigen::SparseMatrix<Eigen::dcomplex> SpMat_cd;

// ===========================================
// type definitions for real OP
   /*typedef double T_scalar;
   typedef Eigen::VectorXd  T_vector;
   typedef Eigen::MatrixXd  T_matrix;
   typedef Eigen::SparseMatrix<double> SpMat_cd;
   */
// ===========================================

// ===========================================
// Boundary condition structure for a single OP
   struct Bound_Cond {
      std::string typeB, typeT, typeL, typeR; // type of BC: D (Dirichlet), N (Neumann), or S (from a 1D Solution)
      double slipB, slipT, slipL, slipR;      // for Neumann BC  ...  the slip length
      double valueB, valueT, valueL, valueR;  // for Dirichlet BC ... the function value

      Bound_Cond& operator= (Bound_Cond& rhs) {
         typeB = rhs.typeB;
         typeT = rhs.typeT;
         typeL = rhs.typeL;
         typeR = rhs.typeR;

         slipB = rhs.slipB;
         slipT = rhs.slipT;
         slipL = rhs.slipL;
         slipR = rhs.slipR;

         valueB = rhs.valueB;
         valueT = rhs.valueT;
         valueL = rhs.valueL;
         valueR = rhs.valueR;
         return *this;
      }
   };

   inline std::ifstream& operator>> (std::ifstream& stream, Bound_Cond& BC) {
      stream >> BC.typeT;  // read type
      stream >> BC.slipT;  // then read both slip length
      stream >> BC.valueT; // and the value

      stream >> BC.typeB;
      stream >> BC.slipB;
      stream >> BC.valueB;

      stream >> BC.typeL;
      stream >> BC.slipL;
      stream >> BC.valueL;

      stream >> BC.typeR;
      stream >> BC.slipR;
      stream >> BC.valueR;

      return stream;
   }

   inline std::ostream& operator<< (std::ostream& out, Bound_Cond& BC) {
      out << "BC.typeT: " << BC.typeT << std::endl;
      if      (BC.typeT == "N") out << "BC.slipT:  " << BC.slipT << std::endl;
      else if (BC.typeT == "D") out << "BC.valueT:  "  << BC.valueT  << std::endl;
      
      out << "BC.typeB: " << BC.typeB << std::endl;
      if      (BC.typeB == "N") out << "BC.slipB:  " << BC.slipB << std::endl;
      else if (BC.typeB == "D") out << "BC.valueB:  "  << BC.valueB  << std::endl;
      
      out << "BC.typeL: " << BC.typeL << std::endl;
      if      (BC.typeL == "N") out << "BC.slipL:  " << BC.slipL << std::endl;
      else if (BC.typeL == "D") out << "BC.valueL:  "  << BC.valueL  << std::endl;
      
      out << "BC.typeR: " << BC.typeR << std::endl;
      if      (BC.typeR == "N") out << "BC.slipR:  " << BC.slipR << std::endl;
      else if (BC.typeR == "D") out << "BC.valueR:  "  << BC.valueR  << std::endl;
      
      return out;
   }
// ===========================================

// ===========================================
// Values used for set up, initial conditions,
// and algorithmic parameters
   struct in_conditions {
	   int Nop;      // number of order parameter components
      int SIZEu;    // the number of points along u-direction (can be x,r, etc)
      int SIZEv;    // ""                         v-direction 
      int SIZEw;    // "" in orthogonal third directions
      double STEP;  // step size
      // physical parameters: temperature and pressure 
      double T; 
      double P; 
      // self-consistency convergence parameters
      int N_loop;   // count limit for iterations
      int maxStore; // max num vectors kept in accerleration
      int wait;     // iterations to wait before using acceleration
      double ACCUR; // the minimum desired accuracy
      double rel_p; // relaxation param

      in_conditions & operator= (in_conditions & rhs) {
         Nop = rhs.Nop;
         SIZEu = rhs.SIZEu;
         SIZEv = rhs.SIZEv;
         SIZEw = rhs.SIZEw;
	      STEP = rhs.STEP;
         T = rhs.T;
         P = rhs.P;
         N_loop = rhs.N_loop;
         maxStore = rhs.maxStore;
         wait = rhs.wait;
         ACCUR = rhs.ACCUR;
         rel_p = rhs.rel_p;
         return *this;
      }
   };
// ===========================================

// ===========================================
// K matrix namespace:
//    to not get the K matrices mixed up with anything else
   namespace kMatrix {
      // a function to give the small Kk matrix to build the solver matrix for smaller systems
      inline Eigen::Matrix2d** smallKMatrix(int Nop, std::vector<bool> D_components) {
         Eigen::Matrix2d** gradK;
         gradK = new Eigen::Matrix2d *[Nop];
         for (int i = 0; i < Nop; i++) gradK[i] = new Eigen::Matrix2d [Nop];

         // make sure the whole starts at 0!
         for (int n = 0; n < Nop; n++)
            for (int m = 0; m < Nop; m++)
               gradK[m][n] << 0, 0, 0, 0;

         for (int n = 0; n < Nop; n++) {
            for (int m = 0; m < Nop; m++) {
               if (D_components[0]) gradK[m][n](0,0) = Kuu[m][n];
               if (D_components[1]) gradK[m][n](0,1) = Kuv[m][n];
               if (D_components[2]) gradK[m][n](1,0) = Kvu[m][n];
               if (D_components[3]) gradK[m][n](1,1) = Kvv[m][n];

               // WILL WE EVER NEED THESE? NO?
               // if (D_components[i]) { // "yy"
               //    //
               // }
               // if (D_components[i]) { // "xy"
               //    //   
               // }
               // if (D_components[i]) { // "yx"
               //    //   
               // }
               // if (D_components[i]) { // "yz"
               //    //
               // }
               // if (D_components[i]) { // "zy"
               //    //
               // }
            }
         }

         return gradK;
      }
   };
// ========================================================

#endif
