#ifndef _complex_includes
#define _complex_includes

// just to keep all the "errors" out while in VS Codes
#include "includes.hpp"

// derived, complex 3-component GL solver class
template<>
class Three_Component_GL_Solver<VectorXcd, dcomplex> : public GL_Solver<VectorXcd, dcomplex>
{
   public:
   Three_Component_GL_Solver(string conditions_file, string boundary_conditions_file)
   {
      ReadConditions(conditions_file);
      ReadBoundaryConditions(boundary_conditions_file);
      size = cond.SIZEu*cond.SIZEv;
   }

   void BuildProblem()
   {
      Build_D_Matrices();
      initialize_OP_matrix();
      SolverMatrix = SolverMatrix_He3Defect();
   }

   void initialize_OP_matrix()
   {
      op_matrix.resize(cond.SIZEv,cond.SIZEu); // initialize matrix
      op_vector.resize(cond.SIZEu*cond.SIZEv*OP_size); // initialize op_vector, for the whole thing (size = num_of_mesh_points * num_OP_components)

      // initialize elements in 'matrix'
      for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            op_matrix(n_v,n_u).initialize(OP_size);
         }
      }

      setVectorForm();// make the vector form available
   }

   // Convert the OP matrix, at all mesh points, into a vector
   void setVectorForm()
   {
      for (int vi = 0; vi < OP_size; vi++) {
         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
               op_vector(ID(cond.SIZEu*cond.SIZEv,n_u,cond.SIZEu,n_v,vi)) = op_matrix(n_v,n_u)(vi);
            }
         }
      }
   }

   // User-defined methods to build the solver matrix and the rhs vector
   SpMat_cd SolverMatrix_He3Defect()
   {
      // to make the code cleaner, define some constants
      double K123 = gl.K1+gl.K2+gl.K3,
               K23  = gl.K2+gl.K3;

      // the matrix to be used by the solver
      SpMat_cd solver_mat(size*OP_size,size*OP_size);

      // initialize each non-zero 'element'
      SpMat_cd elem_00 = Dv2_BD(Axx_BC,0);
      SpMat_cd elem_11 = Dv2_BD(Ayy_BC,1);
      SpMat_cd elem_22 = 3.*Dv2_BD(Azz_BC,2);

      // matrices for placement of each non-zero 'element'
      SpMat_cd M00(OP_size,OP_size), M11(OP_size,OP_size), M22(OP_size,OP_size);
      
      // add a single 1 to each matrix to place the elem_.. matrices in the correct spot
      M00.insert(0,0) = 1.; M11.insert(1,1) = 1.; M22.insert(2,2) = 1.;
      
      // 'place' them using the Kronecker product
      solver_mat = kroneckerProduct( M00, elem_00 ) + kroneckerProduct( M11, elem_11 ) + kroneckerProduct( M22, elem_22 );

      return solver_mat;
   }
   
   VectorXcd RHS_He3Defect()
   {
      int cts = 0;
      Matrix3cd op, op_T, op_dag, op_conj;
      VectorXcd rhs(op_vector.size());
      double Beta_bulk = 6*(gl.B1+gl.B2)+2*(gl.B3+gl.B4+gl.B5);
      Bound_Cond temp_BC;

      // loop through all the OP components in the mesh
      for (int vi = 0; vi < OP_size; vi++) {
         // decide which BC to use
               if (vi == 0) temp_BC = Axx_BC;
         else if (vi == 1) temp_BC = Ayy_BC;
         else if (vi == 2) temp_BC = Azz_BC;
         else cout << "RHS ERROR: OP index out of bounds." << endl;

         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
               int id = ID(cond.SIZEu*cond.SIZEv,n_u,cond.SIZEu,n_v,vi);
               dcomplex val;

               if (temp_BC.typeB == string("Dirichlet") && !n_v) val = temp_BC.bB;
               else if (temp_BC.typeB == string("Neumann") && !n_v) val = 0.;
               else if (temp_BC.typeT == string("Dirichlet") && n_v == cond.SIZEv-1) val = temp_BC.bT;
               else if (temp_BC.typeT == string("Neumann") && n_v == cond.SIZEv-1) val = 0.;
               else
               {
                  auto axx = op_matrix(n_v,n_u)(0);
                  auto azz = op_matrix(n_v,n_u)(2);

                  // assuming all components are real, so a* = a
                  // the Azz_BC: perpendicular
                  if (vi == 2)                 val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*azz + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*azz + 2./5.*pow(abs(azz),2)*azz - azz;
                  // the Axx_BC or Ayy_BC: parallel
                  else if (vi == 0 || vi == 1) val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*axx + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*axx + 2./5.*pow(abs(axx),2)*axx - axx;

                  val *= pow(cond.STEP,2); // because we multiplied the matrix by h^2
               }

               // insert val in rhs, in the matching location
               rhs(id) = val;
            }
         }
      }
      return rhs;
   }

   // given the next guess, update the op at all points on the mesh
   //   and make it available in it's vector form.
   void updateMatrix(VectorXcd& new_guess)
   {
      int sz = cond.SIZEu*cond.SIZEv;
      for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
         for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
            Matrix3cd op;
            op << new_guess(ID(sz,n_u,cond.SIZEu,n_v,0)), 0.,                                     0.,
                  0.,                                     new_guess(ID(sz,n_u,cond.SIZEu,n_v,1)), 0.,
                  0.,                                     0.,                                     new_guess(ID(sz,n_u,cond.SIZEu,n_v,2));
            op_matrix(n_v,n_u).Set_OP(op);
         }
      }
      setVectorForm();
   }

   // use the relaxation method and Anderson Acceleration to solve
   void Solve(VectorXcd& guess)
   {
      cout << "solving..." << endl;
      auto start = std::chrono::system_clock::now();
      VectorXcd f = makeGuess(guess), df(guess.size()); // initialize vectors

      if (no_update.size()) cout << "using no_update" << endl;

      // use LU decomposition to solve the system
      SparseLU<SpMat_cd, COLAMDOrdering<int> > solver;
      solver.analyzePattern(SolverMatrix); // without this, Eigen throws: Eigen::Matrix<int, -1, 1>; ... Assertion `index >= 0 && index < size()' failed.
      solver.factorize(SolverMatrix);

      // check to see if the solver failed
      if (solver.info() == Eigen::Success) cout << "\tSolver: successfully built" << endl;
      else if (solver.info() == Eigen::NumericalIssue) // for debugging non-invertable matrices
      {
         cout << "Solver: numerical issues" << endl;
         // cout << "non-zero indeces of SolverMatrix:\n";
         // for (int k=0; k < SolverMatrix.outerSize(); ++k) for (SpMat_cd::InnerIterator it(SolverMatrix,k); it; ++it) cout << "(" << it.row() << "," << it.col() << ")\t";
         return;
      }

      // time to prepare solving method
      auto end = std::chrono::system_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
      cout << "\ttime: " << elapsed.count() << " seconds." << endl << endl;

      updateMatrix(f); // prepare the matrix for RHS

      int cts = 0; // count loops
      double err;  // to store current error
      VectorXcd rhs = RHS_He3Defect(); // the right hand side

      // the acceleration object
      converg_acceler<VectorXcd> Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update);

      // loop until f converges or until it's gone too long
      do { // use relaxation

         df = solver.solve(rhs)-f; // find the change in f

         if (method == string("acceleration")) Con_Acc.next_vector<MatrixXcd>(f,df,err); // use Anderson Acceleration to converge faster
         else if (method == string("relaxation")) // use normal relaxation
         {
            f += cond.rel_p*df;
            err = df.norm()/f.norm();
         }
         else { cout << "ERROR: Unknown method type given." << endl; return; }

         updateMatrix(f);       // update the matrix for RHS
         rhs = RHS_He3Defect(); // update rhs
         cts++;                 // increment counter

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

   // make an educated guess using the boundary conditions
   VectorXcd makeGuess(VectorXcd& g)
   {
      if (Axx_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1); i < size; i++) g(i) = Axx_BC.bB;
      }
      if (Axx_BC.typeT == string("Dirichlet")) {
         for (int i = 0; i < cond.SIZEu; i++) g(i) = Axx_BC.bT;
      }
      
      if (Ayy_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1)+size; i < 2*size; i++) g(i) = Ayy_BC.bB;
      }
      if (Ayy_BC.typeT == string("Dirichlet")) {
         for (int i = size; i < cond.SIZEu+size; i++) g(i) = Ayy_BC.bT;
      }
      
      if (Azz_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1)+2*size; i < 3*size; i++) g(i) = Azz_BC.bB;
      }
      if (Azz_BC.typeT == string("Dirichlet")) {
         for (int i = 2*size; i < cond.SIZEu+2*size; i++) g(i) = Azz_BC.bT;
      } return g;
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
               if (ls == string("Axx bTop")) {
                  BCs >> Axx_BC.bT;
                  BCs >> Axx_BC.typeT;
               } else if (ls == string("Axx bBott")) {
                  BCs >> Axx_BC.bB;
                  BCs >> Axx_BC.typeB;
               }

               // Ayy
               else if (ls == string("Ayy bTop")) {
                  BCs >> Ayy_BC.bT;
                  BCs >> Ayy_BC.typeT;
               } else if (ls == string("Ayy bBott")) {
                  BCs >> Ayy_BC.bB;
                  BCs >> Ayy_BC.typeB;
               }

               // Azz
               else if (ls == string("Azz bTop")) {
                  BCs >> Azz_BC.bT;
                  BCs >> Azz_BC.typeT;
               } else if (ls == string("Azz bBott")) {
                  BCs >> Azz_BC.bB;
                  BCs >> Azz_BC.typeB;
               }
            }
         }

         BCs.close();
      }
      else cout << "Unable to open file:" << boundary_conditions_file << endl;
   }
   
   // Write all components of the OP, all into one file, of the form:
   //             __x__|__y__|_Axx_|_Axy_| ...
   //              ... | ... | ... | ... | ...
   // separates the components from solution...storing real and imag parts ==> up to 18
   void WriteToFile(VectorXcd& vec, string file_name)
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

   private:
   Matrix<Three_ComponentOrderParam<VectorXcd, dcomplex>,-1,-1> op_matrix;
   VectorXcd op_vector;
};

// derived, complex 3-component GL solver class
template<>
class Five_Component_GL_Solver<VectorXcd, dcomplex> : public GL_Solver<VectorXcd, dcomplex>
{
   public:
   Five_Component_GL_Solver(string conditions_file, string boundary_conditions_file)
   {
      ReadConditions(conditions_file);
      ReadBoundaryConditions(boundary_conditions_file);
      size = cond.SIZEu*cond.SIZEv;
   }

   void BuildProblem()
   {
      Build_D_Matrices();
      initialize_OP_matrix();
      SolverMatrix = SolverMatrix_He3Defect();
   }

   void initialize_OP_matrix()
   {
      op_matrix.resize(cond.SIZEv,cond.SIZEu); // initialize matrix
      op_vector.resize(cond.SIZEu*cond.SIZEv*OP_size); // initialize op_vector, for the whole thing (size = num_of_mesh_points * num_OP_components)

      // initialize elements in 'matrix'
      for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            op_matrix(n_v,n_u).initialize(OP_size);
         }
      }

      setVectorForm();// make the vector form available
   }

   // Convert the OP matrix, at all mesh points, into a vector
   void setVectorForm()
   {
      for (int vi = 0; vi < OP_size; vi++) {
         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
               op_vector(ID(cond.SIZEu*cond.SIZEv,n_u,cond.SIZEu,n_v,vi)) = op_matrix(n_v,n_u)(vi);
            }
         }
      }
   }

   // User-defined methods to build the solver matrix and the rhs vector
   SpMat_cd SolverMatrix_He3Defect()
   {
      // to make the code cleaner, define some constants
      double K123 = gl.K1+gl.K2+gl.K3,
               K23  = gl.K2+gl.K3;

      // the matrix to be used by the solver
      SpMat_cd solver_mat(size*OP_size,size*OP_size);

      // initialize each non-zero 'element'
      SpMat_cd elem_00 = Dv2_BD(Axx_BC,0);
      SpMat_cd elem_11 = Dv2_BD(Ayy_BC,1);
      SpMat_cd elem_22 = 3.*Dv2_BD(Azz_BC,2);

      // matrices for placement of each non-zero 'element'
      SpMat_cd M00(OP_size,OP_size), M11(OP_size,OP_size), M22(OP_size,OP_size);
      
      // add a single 1 to each matrix to place the elem_.. matrices in the correct spot
      M00.insert(0,0) = 1.; M11.insert(1,1) = 1.; M22.insert(2,2) = 1.;
      
      // 'place' them using the Kronecker product
      solver_mat = kroneckerProduct( M00, elem_00 ) + kroneckerProduct( M11, elem_11 ) + kroneckerProduct( M22, elem_22 );

      return solver_mat;
   }
   
   VectorXcd RHS_He3Defect()
   {
      int cts = 0;
      Matrix3cd op, op_T, op_dag, op_conj;
      VectorXcd rhs(op_vector.size());
      double Beta_bulk = 6*(gl.B1+gl.B2)+2*(gl.B3+gl.B4+gl.B5);
      Bound_Cond temp_BC;

      // loop through all the OP components in the mesh
      for (int vi = 0; vi < OP_size; vi++) {
         // decide which BC to use
               if (vi == 0) temp_BC = Axx_BC;
         else if (vi == 1) temp_BC = Ayy_BC;
         else if (vi == 2) temp_BC = Azz_BC;
         else cout << "RHS ERROR: OP index out of bounds." << endl;

         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
               int id = ID(cond.SIZEu*cond.SIZEv,n_u,cond.SIZEu,n_v,vi);
               dcomplex val;

               if (temp_BC.typeB == string("Dirichlet") && !n_v) val = temp_BC.bB;
               else if (temp_BC.typeB == string("Neumann") && !n_v) val = 0.;
               else if (temp_BC.typeT == string("Dirichlet") && n_v == cond.SIZEv-1) val = temp_BC.bT;
               else if (temp_BC.typeT == string("Neumann") && n_v == cond.SIZEv-1) val = 0.;
               else
               {
                  auto axx = op_matrix(n_v,n_u)(0);
                  auto azz = op_matrix(n_v,n_u)(2);

                  // assuming all components are real, so a* = a
                  // the Azz_BC: perpendicular
                  if (vi == 2)                 val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*azz + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*azz + 2./5.*pow(abs(azz),2)*azz - azz;
                  // the Axx_BC or Ayy_BC: parallel
                  else if (vi == 0 || vi == 1) val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*axx + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*axx + 2./5.*pow(abs(axx),2)*axx - axx;

                  val *= pow(cond.STEP,2); // because we multiplied the matrix by h^2
               }

               // insert val in rhs, in the matching location
               rhs(id) = val;
            }
         }
      }
      return rhs;
   }

   // given the next guess, update the op at all points on the mesh
   //   and make it available in it's vector form.
   void updateMatrix(VectorXcd& new_guess)
   {
      int sz = cond.SIZEu*cond.SIZEv;
      for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
         for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
            Matrix3cd op;
            op << new_guess(ID(sz,n_u,cond.SIZEu,n_v,0)), 0.,                                     0.,
                  0.,                                     new_guess(ID(sz,n_u,cond.SIZEu,n_v,1)), 0.,
                  0.,                                     0.,                                     new_guess(ID(sz,n_u,cond.SIZEu,n_v,2));
            op_matrix(n_v,n_u).Set_OP(op);
         }
      }
      setVectorForm();
   }

   // use the relaxation method and Anderson Acceleration to solve
   void Solve(VectorXcd& guess)
   {
      cout << "solving..." << endl;
      auto start = std::chrono::system_clock::now();
      VectorXcd f = makeGuess(guess), df(guess.size()); // initialize vectors

      if (no_update.size()) cout << "using no_update" << endl;

      // use LU decomposition to solve the system
      SparseLU<SpMat_cd, COLAMDOrdering<int> > solver;
      solver.analyzePattern(SolverMatrix); // without this, Eigen throws: Eigen::Matrix<int, -1, 1>; ... Assertion `index >= 0 && index < size()' failed.
      solver.factorize(SolverMatrix);

      // check to see if the solver failed
      if (solver.info() == Eigen::Success) cout << "\tSolver: successfully built" << endl;
      else if (solver.info() == Eigen::NumericalIssue) // for debugging non-invertable matrices
      {
         cout << "Solver: numerical issues" << endl;
         // cout << "non-zero indeces of SolverMatrix:\n";
         // for (int k=0; k < SolverMatrix.outerSize(); ++k) for (SpMat_cd::InnerIterator it(SolverMatrix,k); it; ++it) cout << "(" << it.row() << "," << it.col() << ")\t";
         return;
      }

      // time to prepare solving method
      auto end = std::chrono::system_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
      cout << "\ttime: " << elapsed.count() << " seconds." << endl << endl;

      updateMatrix(f); // prepare the matrix for RHS

      int cts = 0; // count loops
      double err;  // to store current error
      VectorXcd rhs = RHS_He3Defect(); // the right hand side

      // the acceleration object
      converg_acceler<VectorXcd> Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update);

      // loop until f converges or until it's gone too long
      do { // use relaxation

         df = solver.solve(rhs)-f; // find the change in f

         if (method == string("acceleration")) Con_Acc.next_vector<MatrixXcd>(f,df,err); // use Anderson Acceleration to converge faster
         else if (method == string("relaxation")) // use normal relaxation
         {
            f += cond.rel_p*df;
            err = df.norm()/f.norm();
         }
         else { cout << "ERROR: Unknown method type given." << endl; return; }

         updateMatrix(f);       // update the matrix for RHS
         rhs = RHS_He3Defect(); // update rhs
         cts++;                 // increment counter

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

   // make an educated guess using the boundary conditions
   VectorXcd makeGuess(VectorXcd& g)
   {
      if (Axx_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1); i < size; i++) g(i) = Axx_BC.bB;
      }
      if (Axx_BC.typeT == string("Dirichlet")) {
         for (int i = 0; i < cond.SIZEu; i++) g(i) = Axx_BC.bT;
      }
      
      if (Ayy_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1)+size; i < 2*size; i++) g(i) = Ayy_BC.bB;
      }
      if (Ayy_BC.typeT == string("Dirichlet")) {
         for (int i = size; i < cond.SIZEu+size; i++) g(i) = Ayy_BC.bT;
      }
      
      if (Azz_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1)+2*size; i < 3*size; i++) g(i) = Azz_BC.bB;
      }
      if (Azz_BC.typeT == string("Dirichlet")) {
         for (int i = 2*size; i < cond.SIZEu+2*size; i++) g(i) = Azz_BC.bT;
      } return g;
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
               if (ls == string("Axx bTop")) {
                  BCs >> Axx_BC.bT;
                  BCs >> Axx_BC.typeT;
               } else if (ls == string("Axx bBott")) {
                  BCs >> Axx_BC.bB;
                  BCs >> Axx_BC.typeB;
               }

               // Ayy
               else if (ls == string("Ayy bTop")) {
                  BCs >> Ayy_BC.bT;
                  BCs >> Ayy_BC.typeT;
               } else if (ls == string("Ayy bBott")) {
                  BCs >> Ayy_BC.bB;
                  BCs >> Ayy_BC.typeB;
               }

               // Azz
               else if (ls == string("Azz bTop")) {
                  BCs >> Azz_BC.bT;
                  BCs >> Azz_BC.typeT;
               } else if (ls == string("Azz bBott")) {
                  BCs >> Azz_BC.bB;
                  BCs >> Azz_BC.typeB;
               }
            }
         }

         BCs.close();
      }
      else cout << "Unable to open file:" << boundary_conditions_file << endl;
   }
   
   // Write all components of the OP, all into one file, of the form:
   //             __x__|__y__|_Axx_|_Axy_| ...
   //              ... | ... | ... | ... | ...
   // separates the components from solution...storing real and imag parts ==> up to 18
   void WriteToFile(VectorXcd& vec, string file_name)
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

   private:
   Matrix<Five_ComponentOrderParam<VectorXcd, dcomplex>,-1,-1> op_matrix;
   VectorXcd op_vector;
};

#endif