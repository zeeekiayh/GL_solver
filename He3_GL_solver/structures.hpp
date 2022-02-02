#ifndef _structures
#define _structures

#include <iostream>
#include <fstream> // for file in/out

using namespace std;
using std::cout;
using std::endl;
using std::cin;
using std::ostream;
using std::string;

// ===========================================
// Boundary condition structure for a sigle OP
   struct Bound_Cond {
      string typeB, typeT, typeL, typeR; // type of BC: Dirichlet (value) or Neumann (derivative)
      double valB, valT, valL, valR; // for Neumann BC:   slip length
                                     // for Dirichlet BC: function value
      Bound_Cond& operator= (Bound_Cond&);
   };

   ifstream& operator>> (ifstream&, Bound_Cond&);

   ostream& operator<< (ostream&, Bound_Cond&);
// ===========================================

// ===========================================
// Values used for set up, initial conditions,
//   and algorithmic parameters
   struct in_conditions {
      int OP_size;  // number of OP components
      int SIZEu;    // the number of points...
      int SIZEv;    // "" ...
      int SIZEw;    // "" in orthogonal directions
      int N_loop;   // count limit for iterations
      int maxStore; // max num vectors kept in accerleration
      int wait;     // iterations to wait before using acceleration
      double ACCUR; // the minimum desired accuracy
      double STEP;  // step size
      double rel_p; // relaxation param
      double T;     // temperature
      double P;     // pressure
   };
// ===========================================


#endif
