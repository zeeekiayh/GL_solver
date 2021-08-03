#ifndef _real_includes
#define _real_includes

// just to keep all the "errors" out while in VS Codes
#include "includes.hpp"

// derived, real 3-component GL solver class
   template<>
   class Three_Component_GL_Solver<VectorXd, double> : public GL_Solver<VectorXd, double>
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
         SolverMatrix = SolverMatrix_He3Defect(Dv2);
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
      SpMat_d SolverMatrix_He3Defect(SpMat_d& Dv2)
      {
         // to make the code cleaner, define some constants
         double K123 = gl.K1+gl.K2+gl.K3,
                K23  = gl.K2+gl.K3;

         // the matrix to be used by the solver
         SpMat_d solver_mat(Dv2.rows()*OP_size,Dv2.cols()*OP_size);

         // initialize each non-zero 'element'
         SpMat_d elem_00 = Dv2_BD(Axx,0),
                              elem_11 = Dv2_BD(Ayy,1),
                              elem_22 = 3.*Dv2_BD(Azz,2);

         // matrices for placement of each non-zero 'element'
         SpMat_d M00(OP_size,OP_size), M11(OP_size,OP_size), M22(OP_size,OP_size);
         
         // add a single 1 to each matrix to place the elem_.. matrices in the correct spot
         M00.insert(0,0) = 1.; M11.insert(1,1) = 1.; M22.insert(2,2) = 1.;
         
         // 'place' them using the Kronecker product
         solver_mat = kroneckerProduct( M00, elem_00 ) + kroneckerProduct( M11, elem_11 ) + kroneckerProduct( M22, elem_22 );

         return solver_mat;
      }
      
      VectorXd RHS_He3Defect()
      {
         int cts = 0;
         Matrix3d op, op_T, op_dag, op_conj;
         VectorXd rhs(op_vector.size());
         double Beta_bulk = 6*(gl.B1+gl.B2)+2*(gl.B3+gl.B4+gl.B5);
         Bound_Cond temp_BC;

         // loop through all the OP components in the mesh
         for (int vi = 0; vi < OP_size; vi++) {
            // decide which BC to use
                 if (vi == 0) temp_BC = Axx;
            else if (vi == 1) temp_BC = Ayy;
            else if (vi == 2) temp_BC = Azz;
            else cout << "RHS ERROR: OP index out of bounds." << endl;

            for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
               for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
                  int id = ID(cond.SIZEu*cond.SIZEv,n_u,cond.SIZEu,n_v,vi);
                  double val;

                  if (temp_BC.typeB == string("Dirichlet") && !n_v) val = temp_BC.bB;
                  else if (temp_BC.typeB == string("Neumann") && !n_v) val = 0.;
                  else if (temp_BC.typeT == string("Dirichlet") && n_v == cond.SIZEv-1) val = temp_BC.bT;
                  else if (temp_BC.typeT == string("Neumann") && n_v == cond.SIZEv-1) val = 0.;
                  else
                  {
                     auto axx = op_matrix(n_v,n_u)(0);
                     auto azz = op_matrix(n_v,n_u)(2);

                     // assuming all components are real, so a* = a
                     // the Azz: perpendicular
                     if (vi == 2)                 val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*azz + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*azz + 2./5.*pow(abs(azz),2)*azz - azz;
                     // the Axx or Ayy: parallel
                     else if (vi == 0 || vi == 1) val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*axx + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*axx + 2./5.*pow(abs(axx),2)*axx - axx;

                     // auto A = op_matrix(n_v,n_u).GetMatrixForm_He3Defect();
                     // auto A_tran = A.transpose();
                     // auto A_conj = A.conjugate();
                     // auto A_dag = A.adjoint();

                     // val = 2.*gl.B1*(A*A_tran).trace()*A(vi,vi) + 2.*gl.B2*(A*A_dag).trace()*A(vi,vi) + 2.*gl.B3*(A*A_tran*A_conj)(vi,vi) + 2.*gl.B4*(A*A_dag*A)(vi,vi) + 2.*gl.B5*(A_conj*A_tran*A)(vi,vi) + gl.alpha*(A*A_dag).trace();

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
      void updateMatrix(VectorXd& new_guess)
      {
         int sz = cond.SIZEu*cond.SIZEv;
         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
               Matrix3d op;
               op << new_guess(ID(sz,n_u,cond.SIZEu,n_v,0)), 0.,                                     0.,
                     0.,                                     new_guess(ID(sz,n_u,cond.SIZEu,n_v,1)), 0.,
                     0.,                                     0.,                                     new_guess(ID(sz,n_u,cond.SIZEu,n_v,2));
               op_matrix(n_v,n_u).Set_OP(op);
            }
         }
         setVectorForm();
      }

      // use the relaxation method and Anderson Acceleration to solve
      void Solve(VectorXd& guess)
      {
         cout << "solving..." << endl;
         auto start = std::chrono::system_clock::now();
         VectorXd f = makeGuess(guess), df(guess.size()); // initialize vectors

         if (no_update.size()) cout << "using no_update" << endl;

         // use LU decomposition to solve the system
         SparseLU<SpMat_d, COLAMDOrdering<int> > solver;
         solver.analyzePattern(SolverMatrix); // without this, Eigen throws: Eigen::Matrix<int, -1, 1>; ... Assertion `index >= 0 && index < size()' failed.
         solver.factorize(SolverMatrix);

         // check to see if the solver failed
         if (solver.info() == Eigen::Success) cout << "\tSolver: successfully built" << endl;
         else if (solver.info() == Eigen::NumericalIssue) // for debugging non-invertable matrices
         {
            cout << "Solver: numerical issues" << endl;
            // cout << "non-zero indeces of SolverMatrix:\n";
            // for (int k=0; k < SolverMatrix.outerSize(); ++k) for (SpMat_d::InnerIterator it(SolverMatrix,k); it; ++it) cout << "(" << it.row() << "," << it.col() << ")\t";
            return;
         }

         // time to prepare solving method
         auto end = std::chrono::system_clock::now();
         auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
         cout << "\ttime: " << elapsed.count() << " seconds." << endl << endl;

         updateMatrix(f); // prepare the matrix for RHS

         int cts = 0; // count loops
         double err;  // to store current error
         VectorXd rhs = RHS_He3Defect(); // the right hand side

         // the acceleration object
         converg_acceler<VectorXd> Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update);

         // loop until f converges or until it's gone too long
         do { // use relaxation

            df = solver.solve(rhs)-f; // find the change in f

            if (method == string("acceleration")) Con_Acc.next_vector<MatrixXd>(f,df,err); // use Anderson Acceleration to converge faster
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
      VectorXd makeGuess(VectorXd& g)
      {
         if (Axx.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1); i < size; i++) g(i) = Axx.bB;
         }
         if (Axx.typeT == string("Dirichlet")) {
            for (int i = 0; i < cond.SIZEu; i++) g(i) = Axx.bT;
         }
         
         if (Ayy.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1)+size; i < 2*size; i++) g(i) = Ayy.bB;
         }
         if (Ayy.typeT == string("Dirichlet")) {
            for (int i = size; i < cond.SIZEu+size; i++) g(i) = Ayy.bT;
         }
         
         if (Azz.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1)+2*size; i < 3*size; i++) g(i) = Azz.bB;
         }
         if (Azz.typeT == string("Dirichlet")) {
            for (int i = 2*size; i < cond.SIZEu+2*size; i++) g(i) = Azz.bT;
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
            Bound_Cond bc;
            
            while (!BCs.eof()) {
               getline(BCs,line);
               if (line[0] == '#') {           // any good way for error handling here?
                  string ls = line.substr(1);
                  if (ls == string("Axx bTop")) {
                     BCs >> Axx.bT;
                     BCs >> Axx.typeT;
                  } else if (ls == string("Axx bBott")) {
                     BCs >> Axx.bB;
                     BCs >> Axx.typeB;
                  }

                  // Ayy
                  else if (ls == string("Ayy bTop")) {
                     BCs >> Ayy.bT;
                     BCs >> Ayy.typeT;
                  } else if (ls == string("Ayy bBott")) {
                     BCs >> Ayy.bB;
                     BCs >> Ayy.typeB;
                  }

                  // Azz
                  else if (ls == string("Azz bTop")) {
                     BCs >> Azz.bT;
                     BCs >> Azz.typeT;
                  } else if (ls == string("Azz bBott")) {
                     BCs >> Azz.bB;
                     BCs >> Azz.typeB;
                  }
               }
            }

            BCs.close();
         }
         else cout << "Unable to open file:" << boundary_conditions_file << endl;
      }

      double free_energy()
      {
         cout << "free-energy..." << endl;
         if (!solution.size())
         {
            cout << "ERROR: cannot calculate free-energy without a solution." << endl;
            return 0.;
         }

         dcomplex I = 0; // start the integral sum at 0
         dcomplex f_bulk, f_bulk_prev = 0.;
         dcomplex f_grad, f_grad_prev = 0.;
         dcomplex k1, k2, k3; // the 3 derivative values
         Matrix3d A, A_tran, A_dag, A_conj;
         VectorXd I_vals(cond.SIZEv); // value of the integral over distance in z

         Three_ComponentOrderParam<VectorXd, double> A_prev_x, A_next_x, A_prev_z, A_next_z; // the A's for the gradient terms
         cout << "initialized" << endl;

         dcomplex Aajk = 0., Aakj = 0., Aajj = 0., Aakk = 0.;
         // MatrixXd integ(size - 2, size - 2); // value of the integral over the mesh
         Bound_Cond temp_bc_j, temp_bc_k;
         int cts = 0;

         cout << "starting loops..." << endl;

         // loop through all OP's in the mesh
         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {

               // calculate all the needed forms of A
               A = op_matrix(n_v,n_u).GetMatrixForm_He3Defect();
               A_tran = A.transpose();
               A_conj = A.conjugate();
               A_dag = A.adjoint();

               k1 = 0., k2 = 0., k3 = 0.; // start them all at 0 for each OP
               
               if (n_u && n_v && n_u < cond.SIZEu-1 && n_v < cond.SIZEv-1) // if we're not at a boundary,
               {
                  // get all the needed neighbors of A
                  A_prev_x = op_matrix(n_v,n_u-1);
                  A_next_x = op_matrix(n_v,n_u+1);
                  A_prev_z = op_matrix(n_v-1,n_u);
                  A_next_z = op_matrix(n_v+1,n_u);

                  // calculate the gradient at each point:
                  //   loop through the OP at the point
                  for (int a = 0; a < 3; a++) {       // spin index
                     for (int j = 0; j < 3; j++) {    // orbital/derivative index
                        for (int k = 0; k < 3; k++) { // orbital/derivative index

                           // calculate the derivatives depending on the index
                           if (j ==0) { // x-derivative
                              Aajj = (A_next_x(j) - A_prev_x(j))/(cond.STEP*2);
                              Aakj = (A_next_x(k) - A_prev_x(k))/(cond.STEP*2);
                           }
                           // if (j ==1) // y-derivative
                           if (j ==2) { // z-derivative
                              Aajj = (A_next_z(j) - A_prev_z(j))/(cond.STEP*2);
                              Aakj = (A_next_z(k) - A_prev_z(k))/(cond.STEP*2);
                           }

                           if (k == 0) { // x-derivative
                              Aajk = (A_next_x(j) - A_prev_x(j))/(cond.STEP*2);
                              Aakk = (A_next_x(k) - A_prev_x(k))/(cond.STEP*2);
                           }
                           // if (k == 1) // y-derivative
                           if (k == 2) { // z-derivative
                              Aajk = (A_next_z(j) - A_prev_z(j))/(cond.STEP*2);
                              Aakk = (A_next_z(k) - A_prev_z(k))/(cond.STEP*2);
                           }

                           // sum up over the indexes
                           k1 += abs2(Aajk);
                           k2 += conj(Aajj)*Aakk;
                           k3 += conj(Aajk)*Aakj;

                           // reset
                           Aajk = 0., Aakj = 0., Aajj = 0., Aakk = 0.;
                        } // for k
                     } // for j
                  } // for a
               } // if not @ bd
               else { // we're at a boundary

                  // calculate the gradient at each point:
                  //   loop through the OP at the point
                  for (int a = 0; a < 3; a++) {    // spin index
                     for (int j = 0; j < 3; j++) { // orbital/derivative index
                        if (j == 0) temp_bc_j = Axx; // choose BC for the j index
                        if (j == 1) temp_bc_j = Ayy;
                        if (j == 2) temp_bc_j = Azz;

                        for (int k = 0; k < 3; k++) { // orbital/derivative index
                           if (k == 0) temp_bc_k = Axx; // choose BC for the k index
                           if (k == 1) temp_bc_k = Ayy;
                           if (k == 2) temp_bc_k = Azz;

                           // Right now, we only use a step of h, so it is less stable...
                           //    is there a way to work around that?
                           // use BC values to calculate the gradient terms
                           if (j == 0) { // x-derivative
                              if (!n_u &&                temp_bc_j.typeL == string("Neumann")) Aajj = A(a,j)/temp_bc_j.bL;
                              if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Neumann")) Aajj = A(a,j)/temp_bc_j.bR;
                              if (!n_u &&                temp_bc_k.typeL == string("Neumann")) Aakj = A(a,k)/temp_bc_k.bL;
                              if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Neumann")) Aakj = A(a,k)/temp_bc_k.bR;

                              if (!n_u &&                temp_bc_j.typeL == string("Dirichlet")) Aajj = ( op_matrix(n_v,n_u+1)(j) - op_matrix(n_v,n_u)(j) )/(cond.STEP*2);
                              if (!n_u &&                temp_bc_k.typeL == string("Dirichlet")) Aakj = ( op_matrix(n_v,n_u+1)(k) - op_matrix(n_v,n_u)(k) )/(cond.STEP*2);
                              if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Dirichlet")) Aajj = ( op_matrix(n_v,n_u)(j) - op_matrix(n_v,n_u-1)(j) )/(cond.STEP*2);
                              if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Dirichlet")) Aakj = ( op_matrix(n_v,n_u)(k) - op_matrix(n_v,n_u-1)(k) )/(cond.STEP*2);
                           }
                           // if (j == 1) // y derivative
                           if (j == 2) { // z-derivative
                              if (!n_v &&                temp_bc_j.typeB == string("Neumann")) Aajj = A(a,j)/temp_bc_j.bB;
                              if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Neumann")) Aajj = A(a,j)/temp_bc_j.bT;
                              if (!n_v &&                temp_bc_k.typeB == string("Neumann")) Aakj = A(a,k)/temp_bc_k.bB;
                              if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Neumann")) Aakj = A(a,k)/temp_bc_k.bT;

                              if (!n_v &&                temp_bc_j.typeB == string("Dirichlet")) Aajj = ( op_matrix(n_v+1,n_u)(j) - op_matrix(n_v,n_u)(j) )/(cond.STEP*2);
                              if (!n_v &&                temp_bc_k.typeB == string("Dirichlet")) Aakj = ( op_matrix(n_v+1,n_u)(k) - op_matrix(n_v,n_u)(k) )/(cond.STEP*2);
                              if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Dirichlet")) Aajj = ( op_matrix(n_v,n_u)(j) - op_matrix(n_v-1,n_u)(j) )/(cond.STEP*2);
                              if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Dirichlet")) Aakj = ( op_matrix(n_v,n_u)(k) - op_matrix(n_v-1,n_u)(k) )/(cond.STEP*2);
                           }

                           if (k == 0) { // x-derivative
                              if (!n_u &&                temp_bc_j.typeL == string("Neumann")) Aajk = A(a,j)/temp_bc_j.bL;
                              if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Neumann")) Aajk = A(a,j)/temp_bc_j.bR;
                              if (!n_u &&                temp_bc_k.typeL == string("Neumann")) Aakk = A(a,k)/temp_bc_k.bL;
                              if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Neumann")) Aakk = A(a,k)/temp_bc_k.bR;

                              if (!n_u &&                temp_bc_j.typeL == string("Dirichlet")) Aajk = ( op_matrix(n_v,n_u+1)(j) - op_matrix(n_v,n_u)(j) )/(cond.STEP*2);
                              if (!n_u &&                temp_bc_k.typeL == string("Dirichlet")) Aakk = ( op_matrix(n_v,n_u+1)(k) - op_matrix(n_v,n_u)(k) )/(cond.STEP*2);
                              if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Dirichlet")) Aajk = ( op_matrix(n_v,n_u)(j) - op_matrix(n_v,n_u-1)(j) )/(cond.STEP*2);
                              if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Dirichlet")) Aakk = ( op_matrix(n_v,n_u)(k) - op_matrix(n_v,n_u-1)(k) )/(cond.STEP*2);
                           }
                           // if (k == 1) // y derivative
                           if (k == 2) { // z-derivative
                              if (!n_v &&                temp_bc_j.typeB == string("Neumann")) Aajk = A(a,j)/temp_bc_j.bB;
                              if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Neumann")) Aajk = A(a,j)/temp_bc_j.bT;
                              if (!n_v &&                temp_bc_k.typeB == string("Neumann")) Aakk = A(a,k)/temp_bc_k.bB;
                              if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Neumann")) Aakk = A(a,k)/temp_bc_k.bT;

                              if (!n_v &&                temp_bc_j.typeB == string("Dirichlet")) Aajk = ( op_matrix(n_v+1,n_u)(j) - op_matrix(n_v,n_u)(j) )/(cond.STEP*2);
                              if (!n_v &&                temp_bc_k.typeB == string("Dirichlet")) Aakk = ( op_matrix(n_v+1,n_u)(k) - op_matrix(n_v,n_u)(k) )/(cond.STEP*2);
                              if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Dirichlet")) Aajk = ( op_matrix(n_v,n_u)(j) - op_matrix(n_v-1,n_u)(j) )/(cond.STEP*2);
                              if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Dirichlet")) Aakk = ( op_matrix(n_v,n_u)(k) - op_matrix(n_v-1,n_u)(k) )/(cond.STEP*2);
                           }

                           k1 += abs2(Aajk);
                           k2 += conj(Aajj)*Aakk;
                           k3 += conj(Aajk)*Aakj;

                           // reset
                           Aajk = 0., Aakj = 0., Aajj = 0., Aakk = 0.;
                        } // for k
                     } // for j
                  } // for a
               }

               // add them all together to get the gradient term for this OP
               f_grad = gl.K1*k1 + gl.K2*k2 + gl.K3*k3;
               
               f_bulk = gl.B1*abs2((A*A_tran).trace()) + gl.B2*pow((A*A_dag).trace(),2) + gl.B3*(A*A_tran*A_conj*A_dag).trace() + gl.B4*(A*A_dag*A*A_dag).trace() + gl.B5*(A*A_dag*A_conj*A_tran).trace() + gl.alpha*(A*A_dag).trace();
               
               // if (cts) 
               I += ( (f_bulk + f_grad) - 1. )*cond.STEP;

               // f_bulk_prev = f_bulk;
               // f_grad_prev = f_grad;

               cts++;
            }
            I_vals(n_v) = I.real();
         }

         // for (int i = 1; i <= size-2; i++)
         // {
         //    // set the previous values
         //    f_bulk_prev = f_bulk;
         //    f_grad_prev = f_grad;
         //    f_bulk = F_Bulk(i);
         //    // calculate the gradient term
         //    f_grad = F_Grad(i-1,i+1);
         //    // use a rectangular integral approximation, centered at the midpoints
         //    I += ( (f_bulk+f_bulk_prev + f_grad+f_grad_prev)/2. - 1. )*cond.STEP;
         //    integ(i-1) = (f_bulk+f_bulk_prev + f_grad+f_grad_prev)/2.;//I;
         // }

         // save the integrand vector to plot and inspect
         WriteToFile_single_vector(I_vals,"integrand.txt");

         // cout << "The final value of f/f0 = " << integ(integ.size()-1) << endl;
         if (I.imag() >= pow(10,-8)) cout << "WARNING: imaginary part of the free-energy is not zero:" << I.imag() << endl;

         return I.real();
      }


      // Write all components of the OP, all into one file, of the form:
      //             __x__|__y__|_Axx_|_Axy_| ...
      //              ... | ... | ... | ... | ...
      // separates the components from solution...storing real and imag parts ==> up to 18
      void WriteToFile(VectorXd& vec, string file_name)
      {
         std::ofstream data (file_name);
         if (data.is_open())
         {
            for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
               for (int n_u = 0; n_u < cond.SIZEu; n_u++) {

                  string line = to_string(n_u*cond.STEP) + string("\t")
                              + to_string(n_v*cond.STEP) + string("\t"); // add the position components

                  for (int vi = 0; vi < OP_size; vi++) {
                     double element = vec(ID(size,n_u,cond.SIZEu,n_v,vi));
                     line += to_string(element); // add the components of the solution vector
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
               string line = to_string(i*cond.STEP) + string("\t"); // add the position component
               line += to_string(vec(i)); // add the component of the vector
               data << line << endl;
            }
         }
         else cout << "Unable to open file: " << file_name << endl;
      }

      private:
      Matrix<Three_ComponentOrderParam<VectorXd, double>,-1,-1> op_matrix;
      VectorXd op_vector;
   };

#endif