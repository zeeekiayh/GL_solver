#ifndef _complex_includes
#define _complex_includes

// just to keep all the "errors" out while in VS Codes
#include "includes.hpp"

// ==========================================
// Derived class from includes.hpp: GL_Solver
//  for complex, 3-component GL solver
   template<>
   class Three_Component_GL_Solver<dcomplex> : public GL_Solver<dcomplex> {
      public:
      Three_Component_GL_Solver(string conditions_file)
      {
         cout << "here3" << endl;
         ReadConditions(conditions_file);
         cout << "here4" << endl;
         // ReadBoundaryConditions(boundary_conditions_file);
         mSize = cond.SIZEu*cond.SIZEv;
         cout << "here5" << endl;
      }

      // User-defined methods to build the solver matrix and the rhs vector
      SpMat_cd BuildSolverMatrix()
      {
         // the matrix to be used by the solver
         SpMat_cd solver_mat(mSize*OP_size,mSize*OP_size);

         // initialize each non-zero 'element'
         SpMat_cd elem_00 = Dw2_BD(Eta_uu_BC,0);
         SpMat_cd elem_11 = Dw2_BD(Eta_vv_BC,1);
         SpMat_cd elem_22 = 3.*Dw2_BD(Eta_ww_BC,2);

         solver_mat = Place_subMatrix(0,0,OP_size,elem_00) + Place_subMatrix(1,1,OP_size,elem_11) + Place_subMatrix(2,2,OP_size,elem_22);

         return solver_mat;
      }
      
      VectorXcd RHS()
      {
         Matrix3cd op, op_T, op_dag, op_conj;
         VectorXcd rhs(this->op_vector.size());
         double Beta_bulk = 6*(gl.B1+gl.B2)+2*(gl.B3+gl.B4+gl.B5);
         Bound_Cond temp_BC;

         // loop through all the OP components in the mesh
         for (int vi = 0; vi < OP_size; vi++) {
            // decide which BC to use
                 if (vi == 0) temp_BC = Eta_uu_BC;
            else if (vi == 1) temp_BC = Eta_vv_BC;
            else if (vi == 2) temp_BC = Eta_ww_BC;
            else cout << "RHS ERROR: OP index out of bounds." << endl;

            for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
               for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
                  int id = ID(cond.SIZEu*cond.SIZEv,n_u,cond.SIZEu,n_v,vi);
                  dcomplex val;

                       if (temp_BC.typeB == string("Dirichlet") && n_v == 0) val = temp_BC.valB;
                  else if (temp_BC.typeT == string("Dirichlet") && n_v == cond.SIZEv-1) val = temp_BC.valT;
                  else
                  {
                     auto a = matrix_operator(op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size);
                     auto axx = a(0);
                     auto azz = a(2);

                     // the Eta_ww_BC: perpendicular
                     if (vi == 2)                 val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*conj(azz) + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*azz + 2./5.*pow(abs(azz),2)*azz - azz;
                     // the Eta_uu_BC or Eta_vv_BC: parallel
                     else if (vi == 0 || vi == 1) val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*conj(axx) + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*axx + 2./5.*pow(abs(axx),2)*axx - axx;


                     val *= pow(cond.STEP,2); // because we multiplied the matrix by h^2
                  }

                  // insert val in rhs, in the matching location
                  rhs(id) = val;
               }
            }
         }
         return rhs;
      }

      // make an educated guess using the boundary conditions
      VectorXcd makeGuess(VectorXcd& g)
      {
         // Eta_uu (xx) BC's
         if (Eta_uu_BC.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1); i < mSize; i++) g(i) = Eta_uu_BC.valB;
         }
         if (Eta_uu_BC.typeT == string("Dirichlet")) {
            for (int i = 0; i < cond.SIZEu; i++) g(i) = Eta_uu_BC.valT;
         }
         
         // Eta_vv (yy) BC's
         if (Eta_vv_BC.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1)+mSize; i < 2*mSize; i++) g(i) = Eta_vv_BC.valB;
         }
         if (Eta_vv_BC.typeT == string("Dirichlet")) {
            for (int i = mSize; i < cond.SIZEu+mSize; i++) g(i) = Eta_vv_BC.valT;
         }
         
         // Eta_ww (zz) BC's
         if (Eta_ww_BC.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1)+2*mSize; i < 3*mSize; i++) g(i) = Eta_ww_BC.valB;
         }
         if (Eta_ww_BC.typeT == string("Dirichlet")) {
            for (int i = 2*mSize; i < cond.SIZEu+2*mSize; i++) g(i) = Eta_ww_BC.valT;
         } return g;
      }

      double free_energy()
      {
         if (!solution.size())
         {
            cout << "ERROR: cannot calculate free-energy without a solution." << endl;
            return 0.;
         }

         double beta_B = 6.*(gl.B1+gl.B2) + 2.*(gl.B3+gl.B4+gl.B5);
         dcomplex I = 0.; // start the integral sum at 0
         dcomplex f_bulk = 0.;
         dcomplex f_grad = 0.;
         dcomplex k1, k2, k3; // the 3 derivative values
         Matrix3cd Eta_, Eta__tran, Eta__dag, Eta__conj;
         VectorXd I_vals(cond.SIZEv); // value of the integral over distance in z
         VectorXcd I_bulk(mSize), I_grad(mSize);
         OrderParam<dcomplex> Eta__prev_x, Eta__next_x, Eta__prev_z, Eta__next_z; // the Eta_'s for the gradient terms
         dcomplex Eta_ajk = 0., Eta_akj = 0., Eta_ajj = 0., Eta_akk = 0.;
         Bound_Cond temp_bc_j, temp_bc_k;

         // loop through all OP's in the mesh
         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {

               // calculate all the needed forms of Eta_
               Eta_ = matrix_operator(op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size).GetMatrixForm();
               Eta__tran = Eta_.transpose();
               Eta__conj = Eta_.conjugate();
               Eta__dag = Eta_.adjoint();

               k1 = 0., k2 = 0., k3 = 0.; // start them all at 0 for each OP
               
               if (n_u && n_v && n_u < cond.SIZEu-1 && n_v < cond.SIZEv-1) // if we're not at a boundary,
               {
                  // get all the needed neighbors of Eta_
                  Eta__prev_x = matrix_operator(op_vector,n_v,n_u-1,cond.SIZEv,cond.SIZEu,OP_size);
                  Eta__next_x = matrix_operator(op_vector,n_v,n_u+1,cond.SIZEv,cond.SIZEu,OP_size);
                  Eta__prev_z = matrix_operator(op_vector,n_v-1,n_u,cond.SIZEv,cond.SIZEu,OP_size);
                  Eta__next_z = matrix_operator(op_vector,n_v+1,n_u,cond.SIZEv,cond.SIZEu,OP_size);

                  // calculate the gradient at each point:
                  //   loop through the OP at the point
                  for (int a = 0; a < 3; a++) {       // spin index
                     for (int j = 0; j < 3; j++) {    // orbital/derivative index
                        for (int k = 0; k < 3; k++) { // orbital/derivative index

                           // calculate the derivatives depending on the index.
                           // divide by h^2 later for nicer code.
                           if (j ==0) { // x-derivative
                              Eta_ajj = (Eta__next_x(j) - Eta__prev_x(j))/cond.STEP*2.;
                              Eta_akj = (Eta__next_x(k) - Eta__prev_x(k))/cond.STEP*2.;
                           }
                           // if (j ==1) // y-derivative
                           if (j ==2) { // z-derivative
                              Eta_ajj = (Eta__next_z(j) - Eta__prev_z(j))/cond.STEP*2.;
                              Eta_akj = (Eta__next_z(k) - Eta__prev_z(k))/cond.STEP*2.;
                           }

                           if (k == 0) { // x-derivative
                              Eta_ajk = (Eta__next_x(j) - Eta__prev_x(j))/cond.STEP*2.;
                              Eta_akk = (Eta__next_x(k) - Eta__prev_x(k))/cond.STEP*2.;
                           }
                           // if (k == 1) // y-derivative
                           if (k == 2) { // z-derivative
                              Eta_ajk = (Eta__next_z(j) - Eta__prev_z(j))/cond.STEP*2.;
                              Eta_akk = (Eta__next_z(k) - Eta__prev_z(k))/cond.STEP*2.;
                           }

                           // sum up over the indexes
                           k1 += norm(Eta_ajk);
                           k2 += conj(Eta_ajj)*Eta_akk;
                           k3 += conj(Eta_ajk)*Eta_akj;

                           // reset
                           Eta_ajk = 0., Eta_akj = 0., Eta_ajj = 0., Eta_akk = 0.;
                        } // for k
                     } // for j
                  } // for a
               } // if not @ bd
               else { // we're at a boundary
                  OrderParam<dcomplex> Eta__vu, Eta__vup, Eta__vum, Eta__vpu, Eta__vmu;

                  // calculate the gradient at each point:
                  //   loop through the OP at the point
                  for (int a = 0; a < 3; a++) {    // spin index
                     for (int j = 0; j < 3; j++) { // orbital/derivative index
                        if (j == 0) temp_bc_j = Eta_uu_BC; // choose BC for the j index
                        if (j == 1) temp_bc_j = Eta_ww_BC;
                        if (j == 2) temp_bc_j = Eta_vv_BC;

                        for (int k = 0; k < 3; k++) { // orbital/derivative index
                           if (k == 0) temp_bc_k = Eta_uu_BC; // choose BC for the k index
                           if (k == 1) temp_bc_k = Eta_ww_BC;
                           if (k == 2) temp_bc_k = Eta_vv_BC;

                           Eta__vu = matrix_operator(op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size);
                           if (n_u != cond.SIZEu-1) Eta__vup = matrix_operator(op_vector,n_v,n_u+1,cond.SIZEv,cond.SIZEu,OP_size); // 'up' = u plus = u+1
                           if (n_u != 0) Eta__vum = matrix_operator(op_vector,n_v,n_u-1,cond.SIZEv,cond.SIZEu,OP_size); // 'um' = u minus = u-1
                           if (n_v != cond.SIZEv-1) Eta__vpu = matrix_operator(op_vector,n_v+1,n_u,cond.SIZEv,cond.SIZEu,OP_size); // 'vp' = v plus = v+1
                           if (n_v != 0) Eta__vmu = matrix_operator(op_vector,n_v-1,n_u,cond.SIZEv,cond.SIZEu,OP_size); // 'vm' = v minus = v-1

                           // Right now, we only use a step of h, so it is less stable...
                           //    is there a way to work around that?
                           // use BC values to calculate the gradient terms
                           if (j == 0) { // x-derivative
                              if (n_u == 0 &&            temp_bc_j.typeL == string("Neumann")) Eta_ajj = Eta_(a,j)/temp_bc_j.valL;
                              if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Neumann")) Eta_ajj = Eta_(a,j)/temp_bc_j.valR;
                              if (n_u == 0 &&            temp_bc_k.typeL == string("Neumann")) Eta_akj = Eta_(a,k)/temp_bc_k.valL;
                              if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Neumann")) Eta_akj = Eta_(a,k)/temp_bc_k.valR;

                              if (n_u == 0 &&            temp_bc_j.typeL == string("Dirichlet")) Eta_ajj = ( Eta__vup(j) - Eta__vu(j) )/cond.STEP;
                              if (n_u == 0 &&            temp_bc_k.typeL == string("Dirichlet")) Eta_akj = ( Eta__vup(k) - Eta__vu(k) )/cond.STEP;
                              if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Dirichlet")) Eta_ajj = ( Eta__vu(j) - Eta__vum(j) )/cond.STEP;
                              if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Dirichlet")) Eta_akj = ( Eta__vu(k) - Eta__vum(k) )/cond.STEP;
                           }
                           // if (j == 1) // y derivative
                           if (j == 2) { // z-derivative
                              if (n_v == 0 &&            temp_bc_j.typeB == string("Neumann")) Eta_ajj = Eta_(a,j)/temp_bc_j.valB;
                              if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Neumann")) Eta_ajj = Eta_(a,j)/temp_bc_j.valT;
                              if (n_v == 0 &&            temp_bc_k.typeB == string("Neumann")) Eta_akj = Eta_(a,k)/temp_bc_k.valB;
                              if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Neumann")) Eta_akj = Eta_(a,k)/temp_bc_k.valT;

                              if (n_v == 0 &&            temp_bc_j.typeB == string("Dirichlet")) Eta_ajj = ( Eta__vpu(j) - Eta__vu(j) )/cond.STEP;
                              if (n_v == 0 &&            temp_bc_k.typeB == string("Dirichlet")) Eta_akj = ( Eta__vpu(k) - Eta__vu(k) )/cond.STEP;
                              if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Dirichlet")) Eta_ajj = ( Eta__vu(j) - Eta__vmu(j) )/cond.STEP;
                              if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Dirichlet")) Eta_akj = ( Eta__vu(k) - Eta__vmu(k) )/cond.STEP;
                           }

                           if (k == 0) { // x-derivative
                              if (n_u == 0 &&            temp_bc_j.typeL == string("Neumann")) Eta_ajk = Eta_(a,j)/temp_bc_j.valL;
                              if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Neumann")) Eta_ajk = Eta_(a,j)/temp_bc_j.valR;
                              if (n_u == 0 &&            temp_bc_k.typeL == string("Neumann")) Eta_akk = Eta_(a,k)/temp_bc_k.valL;
                              if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Neumann")) Eta_akk = Eta_(a,k)/temp_bc_k.valR;

                              if (n_u == 0 &&            temp_bc_j.typeL == string("Dirichlet")) Eta_ajk = ( Eta__vup(j) - Eta__vu(j) )/cond.STEP;
                              if (n_u == 0 &&            temp_bc_k.typeL == string("Dirichlet")) Eta_akk = ( Eta__vup(k) - Eta__vu(k) )/cond.STEP;
                              if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Dirichlet")) Eta_ajk = ( Eta__vu(j) - Eta__vum(j) )/cond.STEP;
                              if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Dirichlet")) Eta_akk = ( Eta__vu(k) - Eta__vum(k) )/cond.STEP;
                           }
                           // if (k == 1) // y derivative
                           if (k == 2) { // z-derivative
                              if (n_v == 0 &&            temp_bc_j.typeB == string("Neumann")) Eta_ajk = Eta_(a,j)/temp_bc_j.valB;
                              if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Neumann")) Eta_ajk = Eta_(a,j)/temp_bc_j.valT;
                              if (n_v == 0 &&            temp_bc_k.typeB == string("Neumann")) Eta_akk = Eta_(a,k)/temp_bc_k.valB;
                              if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Neumann")) Eta_akk = Eta_(a,k)/temp_bc_k.valT;

                              if (n_v == 0 &&            temp_bc_j.typeB == string("Dirichlet")) Eta_ajk = ( Eta__vpu(j) - Eta__vu(j) )/cond.STEP;
                              if (n_v == 0 &&            temp_bc_k.typeB == string("Dirichlet")) Eta_akk = ( Eta__vpu(k) - Eta__vu(k) )/cond.STEP;
                              if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Dirichlet")) Eta_ajk = ( Eta__vu(j) - Eta__vmu(j) )/cond.STEP;
                              if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Dirichlet")) Eta_akk = ( Eta__vu(k) - Eta__vmu(k) )/cond.STEP;
                           }

                           k1 += norm(Eta_ajk);
                           k2 += conj(Eta_ajj)*Eta_akk;
                           k3 += conj(Eta_ajk)*Eta_akj;

                           // reset
                           Eta_ajk = 0., Eta_akj = 0., Eta_ajj = 0., Eta_akk = 0.;
                        } // for k
                     } // for j
                  } // for a
               }

               // add them all together to get the gradient term for this OP
               f_grad = -2./3.*(k1 + k2 + k3);
               // f_grad = gl.K1*k1 + gl.K2*k2 + gl.K3*k3; // the general way
               I_grad(ID(mSize,n_u,cond.SIZEu,n_v,0)) = f_grad;
               
               // this value is not normalized, while the values put into here have been...this can only be used if the GL equations we solve are not normalized
               // f_bulk = gl.B1*norm((Eta_*Eta__tran).trace()) + gl.B2*pow((Eta_*Eta__dag).trace(),2) + gl.B3*(Eta_*Eta__tran*Eta__conj*Eta__dag).trace() + gl.B4*(Eta_*Eta__dag*Eta_*Eta__dag).trace() + gl.B5*(Eta_*Eta__dag*Eta__conj*Eta__tran).trace() + gl.alpha*(Eta_*Eta__dag).trace();
               
               f_bulk = (gl.B1*norm((Eta_ * Eta__tran).trace())
                        +gl.B2*pow( (Eta_ * Eta__dag).trace(), 2)
                        +gl.B3*(Eta_ * Eta__tran * Eta__conj * Eta__dag).trace()
                        +gl.B4*(Eta_ * Eta__dag * Eta_ * Eta__dag).trace()
                        +gl.B5*(Eta_ * Eta__dag * Eta__conj * Eta__tran).trace()
                        )*-1./(beta_B*9) + 2./3.*(Eta_ * Eta__dag).trace();
               I_bulk(ID(mSize,n_u,cond.SIZEu,n_v,0)) = f_bulk;

               I += ( (f_bulk + f_grad) - 1. )*cond.STEP;
            }
            I_vals(n_v) = I.real();
         }

         // save the integrand vector to plot and inspect
         WriteToFile_single_vector(I_vals,"integrand.txt");
         WriteToFile_single_vector(I_bulk,"integ_bulk.txt");
         WriteToFile_single_vector(I_grad,"integ_grad.txt");

         //cout << "The final value of f/f0 = " << integ(integ.size()-1) << endl;
         if (I.imag() >= pow(10,-8)) cout << "WEta_RNING: imaginary part of the free-energy is not zero:" << I.imag() << endl;

         return I.real();
      }
   };
// ==========================================


// ==========================================
// Derived class from includes.hpp: GL_Solver
//  for complex, 5-component GL solver
   template<>
   class Five_Component_GL_Solver<dcomplex> : public GL_Solver<dcomplex> {
      public:
      Five_Component_GL_Solver(string conditions_file)
      {
         ReadConditions(conditions_file);
         // ReadBoundaryConditions(boundary_conditions_file);
         mSize = cond.SIZEu*cond.SIZEv;
      }

      // User-defined methods to build the solver matrix and the rhs vector
      SpMat_cd BuildSolverMatrix()
      {
         // cout << "\tin BuildSolverMatrix()" << endl;
         // to make the code cleaner, define some constants
         double K23  = gl.K2+gl.K3;
         double K123 = gl.K1+gl.K2+gl.K3;

         // the matrix to be used by the solver
         SpMat_cd solver_mat(mSize*OP_size,mSize*OP_size);
         // cout << "\t\tsolver_mat initialized" << endl;

         // initialize each non-zero 'element'
         SpMat_cd elem_00 = K123*Du2_BD(Eta_uu_BC,0) + gl.K1*Dw2_BD(Eta_uu_BC,0);
         SpMat_cd elem_01 = K23*Duw_BD(Eta_uw_BC,0);
         SpMat_cd elem_10 = K23*Duw_BD(Eta_uu_BC,1);
         SpMat_cd elem_11 = K123*Dw2_BD(Eta_uw_BC,1) + gl.K1*Du2_BD(Eta_uw_BC,1);
         SpMat_cd elem_22 = gl.K1*(Du2_BD(Eta_vv_BC,2) + Dw2_BD(Eta_vv_BC,2));
         SpMat_cd elem_33 = K123*Du2_BD(Eta_wu_BC,3) + gl.K1*Dw2_BD(Eta_wu_BC,3);
         SpMat_cd elem_34 = K23*Duw_BD(Eta_ww_BC,3);
         SpMat_cd elem_43 = K23*Duw_BD(Eta_wu_BC,4);
         SpMat_cd elem_44 = K123*Dw2_BD(Eta_ww_BC,4) + gl.K1*Du2_BD(Eta_ww_BC,4);
         // cout << "\t\t'elements' initialized" << endl;

         // cout << "\t\tusing 'Place_subMatrix()'" << endl;
         // insert each 'element' from above into their respective locations
         //    and add them together to get the solver matrix
         solver_mat = Place_subMatrix(0,0,OP_size,elem_00) + Place_subMatrix(0,1,OP_size,elem_01)
                     + Place_subMatrix(1,0,OP_size,elem_10) + Place_subMatrix(1,1,OP_size,elem_11) 
                     + Place_subMatrix(2,2,OP_size,elem_22) + Place_subMatrix(3,3,OP_size,elem_33)
                     + Place_subMatrix(3,4,OP_size,elem_34) + Place_subMatrix(4,3,OP_size,elem_43)
                     + Place_subMatrix(4,4,OP_size,elem_44);
         // cout << "\t\tusing 'Place_subMatrix()': success" << endl;

         return solver_mat;
      }
      
      VectorXcd RHS()
      {
         Matrix3cd op, op_T, op_dag, op_conj;
         VectorXcd rhs(this->op_vector.size());
         Bound_Cond temp_BC;

         // define common values for brevity
         double beta_bulk = 5.*(gl.B1+gl.B2)+17./5.*(gl.B3+gl.B4+gl.B5);
         double beta_1_5 = gl.B1 + gl.B2 + gl.B3 + gl.B4 + gl.B5;
         double beta_234 = gl.B2 + gl.B3 + gl.B4;
         double beta_245 = gl.B2 + gl.B4 + gl.B5;
         double beta_13 = gl.B1 + gl.B3;
         double beta_15 = gl.B1 + gl.B5;

         // loop through all the OP components in the mesh
         for (int vi = 0; vi < OP_size; vi++) {
            // decide which BC to use
                 if (vi == 0) temp_BC = Eta_uu_BC;
            else if (vi == 1) temp_BC = Eta_uw_BC;
            else if (vi == 2) temp_BC = Eta_vv_BC;
            else if (vi == 3) temp_BC = Eta_wu_BC;
            else if (vi == 4) temp_BC = Eta_ww_BC;
            else cout << "RHS ERROR: OP index out of bounds." << endl;

            for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
               for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
                  int id = ID(cond.SIZEu*cond.SIZEv,n_u,cond.SIZEu,n_v,vi);
                  dcomplex val;

                       if (temp_BC.typeB == string("Dirichlet") && n_v == 0) val = temp_BC.valB;
                  else if (temp_BC.typeT == string("Dirichlet") && n_v == cond.SIZEv-1) val = temp_BC.valT;
                  else
                  {
                     auto a = matrix_operator(this->op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size);
                     auto axx = a(0);
                     auto axz = a(1);
                     auto ayy = a(2);
                     auto azx = a(3);
                     auto azz = a(4);

                     switch (vi)
                     {
                     case 0: // Eta_xx                                                                   b2 here?            and                b3 here?
                        val = axx + axx/beta_bulk * 2.*( beta_1_5*norm(axx) + beta_245*norm(axz) + gl.B2*norm(ayy) + beta_234*norm(azx) + gl.B3*norm(azz) )
                              + 2./beta_bulk*( gl.B4*axz*azx*conj(azz) + gl.B3*axz*conj(azx)*azz + gl.B5*conj(axz)*azx*azz )
                              + 2.*conj(axx)/beta_bulk*( beta_13*pow(axz,2) + gl.B1*(pow(ayy,2)+pow(azz,2)) + beta_15*pow(azx,2) );
                        break;
                     case 1: // Eta_xz
                        val = axz + axz/beta_bulk * 2.*( beta_245*norm(axx) + beta_1_5*norm(axz) + gl.B2*norm(ayy) + gl.B2*norm(azx) + beta_234*norm(azz) )
                              + 2./beta_bulk*( gl.B3*axx*azx*conj(azz) + gl.B4*axx*conj(azx)*azz + gl.B5*conj(axx)*azx*azz )
                              + 2.*conj(axz)/beta_bulk*( beta_13*pow(axx,2) + gl.B1*(pow(ayy,2)+pow(azx,2)) + beta_15*pow(azz,2) );
                        break;
                     case 2: // Eta_yy
                        val = ayy + 2.*ayy*gl.B2/beta_bulk * (norm(axx) + norm(axz) + norm(azx) + norm(azz)) + 2.*beta_1_5/beta_bulk*ayy*norm(ayy)
                              + 2.*gl.B1*conj(ayy)/beta_bulk*(pow(axx,2) + pow(axz,2) + pow(azx,2) + pow(azz,2));
                        break;
                     case 3: // Eta_zx
                        val = azx + azx/beta_bulk * 2.*( beta_234*norm(axx) + gl.B2*norm(axz) + gl.B2*norm(ayy) + beta_1_5*norm(azx) + beta_245*norm(azz) )
                              + 2./beta_bulk*( gl.B5*axx*axz*conj(azz) + gl.B4*axx*conj(axz)*azz + gl.B3*conj(axx)*axz*azz )
                              + 2.*conj(azx)/beta_bulk*( beta_15*pow(axx,2) + gl.B1*(pow(axz,2)+pow(ayy,2)) + beta_13*pow(azx,2) );
                        break;
                     case 4: // Eta_zz
                        val = azz + azz/beta_bulk * 2.*( gl.B2*norm(axx) + beta_234*norm(axz) + gl.B2*norm(ayy) + beta_245*norm(azx) + beta_1_5*norm(azz) )
                              + 2./beta_bulk*( gl.B5*axx*azx*conj(azx) + gl.B3*axx*conj(azx)*azx + gl.B4*conj(axx)*azx*azx )
                              + 2.*conj(azz)/beta_bulk*( beta_15*pow(axz,2) + gl.B1*(pow(axx,2)+pow(ayy,2)) + beta_13*pow(azx,2) );
                        break;
                     
                     default:
                        break;
                     }

                     val *= pow(cond.STEP,2); // because we multiplied the matrix by h^2
                  }

                  // insert val in rhs, in the matching location
                  rhs(id) = val;
               }
            }
         }
         return rhs;
      }

      // make an educated guess using the boundary conditions
      VectorXcd makeGuess(VectorXcd& g)
      {
         if (Eta_uu_BC.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1); i < mSize; i++) g(i) = Eta_uu_BC.valB;
         }
         if (Eta_uu_BC.typeT == string("Dirichlet")) {
            for (int i = 0; i < cond.SIZEu; i++) g(i) = Eta_uu_BC.valT;
         }
         
         if (Eta_vv_BC.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1)+mSize; i < 2*mSize; i++) g(i) = Eta_vv_BC.valB;
         }
         if (Eta_vv_BC.typeT == string("Dirichlet")) {
            for (int i = mSize; i < cond.SIZEu+mSize; i++) g(i) = Eta_vv_BC.valT;
         }
         
         if (Eta_ww_BC.typeB == string("Dirichlet")) {
            for (int i = cond.SIZEu*(cond.SIZEv-1)+2*mSize; i < 3*mSize; i++) g(i) = Eta_ww_BC.valB;
         }
         if (Eta_ww_BC.typeT == string("Dirichlet")) {
            for (int i = 2*mSize; i < cond.SIZEu+2*mSize; i++) g(i) = Eta_ww_BC.valT;
         } return g;
      }
   };
// ==========================================

#endif