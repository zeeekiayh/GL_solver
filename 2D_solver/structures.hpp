#ifndef structures_hpp_
#define structures_hpp_

#include <iostream>
#include <fstream> // for file in/out
#include <complex> 
#include <typeinfo>

#include <eigen/Eigen/Dense>
#include <eigen/Eigen/Sparse>

// ===========================================
// type definitions for complex OP

typedef std::complex<double> T_scalar;
typedef Eigen::VectorXcd  T_vector;
typedef Eigen::MatrixXcd  T_matrix;
typedef Eigen::SparseMatrix<Eigen::dcomplex> SpMat_cd;

// ===========================================
// type definitions for real OP

/*
typedef double T_scalar;
typedef Eigen::VectorXd  T_vector;
typedef Eigen::MatrixXd  T_matrix;
typedef Eigen::SparseMatrix<double> SpMat_cd;
*/

// ===========================================
// Boundary condition structure for a single OP
   struct Bound_Cond {
      std::string typeB, typeT, typeL, typeR; // type of BC: Dirichlet (value) or Neumann (derivative)
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


#endif
