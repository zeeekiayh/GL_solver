#ifndef _complex_includes
#define _complex_includes

#include "includes.hpp"

// derived, complex 3-component GL solver class
template<>
class Three_Component_GL_Solver<dcomplex> : public GL_Solver<dcomplex>
{
   public:
   Three_Component_GL_Solver(string conditions_file, string boundary_conditions_file)
   {
      ReadConditions(conditions_file);
      ReadBoundaryConditions(boundary_conditions_file);
      size = cond.SIZEu*cond.SIZEv;
   }

   // User-defined methods to build the solver matrix and the rhs vector
   SpMat_cd BuildSolverMatrix()
   {
      // the matrix to be used by the solver
      SpMat_cd solver_mat(size*OP_size,size*OP_size);

      // initialize each non-zero 'element'
      SpMat_cd elem_00 = Dw2_BD(Auu_BC,0);
      SpMat_cd elem_11 = Dw2_BD(Avv_BC,1);
      SpMat_cd elem_22 = 3.*Dw2_BD(Aww_BC,2);

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
              if (vi == 0) temp_BC = Auu_BC;
         else if (vi == 1) temp_BC = Avv_BC;
         else if (vi == 2) temp_BC = Aww_BC;
         else cout << "RHS ERROR: OP index out of bounds." << endl;

         for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
            for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
               int id = ID(cond.SIZEu*cond.SIZEv,n_u,cond.SIZEu,n_v,vi);
               dcomplex val;

                    if (temp_BC.typeB == string("Dirichlet") && n_v == 0) val = temp_BC.bB;
               else if (temp_BC.typeT == string("Dirichlet") && n_v == cond.SIZEv-1) val = temp_BC.bT;
               else
               {
                  auto a = matrix_operator(op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size);
                  auto axx = a(0);
                  auto azz = a(2);

                  // the Aww_BC: perpendicular
                  if (vi == 2)                 val = -1./5.*( 2.*pow(axx,2)+pow(azz,2) )*conj(azz) + 2./5.*( 2.*pow(abs(axx),2)+pow(abs(azz),2) )*azz + 2./5.*pow(abs(azz),2)*azz - azz;
                  // the Auu_BC or Avv_BC: parallel
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
      // Auu (xx) BC's
      if (Auu_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1); i < size; i++) g(i) = Auu_BC.bB;
      }
      if (Auu_BC.typeT == string("Dirichlet")) {
         for (int i = 0; i < cond.SIZEu; i++) g(i) = Auu_BC.bT;
      }
      
      // Avv (yy) BC's
      if (Avv_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1)+size; i < 2*size; i++) g(i) = Avv_BC.bB;
      }
      if (Avv_BC.typeT == string("Dirichlet")) {
         for (int i = size; i < cond.SIZEu+size; i++) g(i) = Avv_BC.bT;
      }
      
      // Aww (zz) BC's
      if (Aww_BC.typeB == string("Dirichlet")) {
         for (int i = cond.SIZEu*(cond.SIZEv-1)+2*size; i < 3*size; i++) g(i) = Aww_BC.bB;
      }
      if (Aww_BC.typeT == string("Dirichlet")) {
         for (int i = 2*size; i < cond.SIZEu+2*size; i++) g(i) = Aww_BC.bT;
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
      Matrix3cd A, A_tran, A_dag, A_conj;
      VectorXd I_vals(cond.SIZEv); // value of the integral over distance in z
      VectorXcd I_bulk(size), I_grad(size);
      OrderParam<dcomplex> A_prev_x, A_next_x, A_prev_z, A_next_z; // the A's for the gradient terms
      dcomplex Aajk = 0., Aakj = 0., Aajj = 0., Aakk = 0.;
      Bound_Cond temp_bc_j, temp_bc_k;

      // loop through all OP's in the mesh
      for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
         for (int n_u = 0; n_u < cond.SIZEu; n_u++) {

            // calculate all the needed forms of A
            A = matrix_operator(op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size).GetMatrixForm();
            A_tran = A.transpose();
            A_conj = A.conjugate();
            A_dag = A.adjoint();

            k1 = 0., k2 = 0., k3 = 0.; // start them all at 0 for each OP
            
            if (n_u && n_v && n_u < cond.SIZEu-1 && n_v < cond.SIZEv-1) // if we're not at a boundary,
            {
               // get all the needed neighbors of A
               A_prev_x = matrix_operator(op_vector,n_v,n_u-1,cond.SIZEv,cond.SIZEu,OP_size);
               A_next_x = matrix_operator(op_vector,n_v,n_u+1,cond.SIZEv,cond.SIZEu,OP_size);
               A_prev_z = matrix_operator(op_vector,n_v-1,n_u,cond.SIZEv,cond.SIZEu,OP_size);
               A_next_z = matrix_operator(op_vector,n_v+1,n_u,cond.SIZEv,cond.SIZEu,OP_size);

               // calculate the gradient at each point:
               //   loop through the OP at the point
               for (int a = 0; a < 3; a++) {       // spin index
                  for (int j = 0; j < 3; j++) {    // orbital/derivative index
                     for (int k = 0; k < 3; k++) { // orbital/derivative index

                        // calculate the derivatives depending on the index.
                        // divide by h^2 later for nicer code.
                        if (j ==0) { // x-derivative
                           Aajj = (A_next_x(j) - A_prev_x(j))/cond.STEP*2.;
                           Aakj = (A_next_x(k) - A_prev_x(k))/cond.STEP*2.;
                        }
                        // if (j ==1) // y-derivative
                        if (j ==2) { // z-derivative
                           Aajj = (A_next_z(j) - A_prev_z(j))/cond.STEP*2.;
                           Aakj = (A_next_z(k) - A_prev_z(k))/cond.STEP*2.;
                        }

                        if (k == 0) { // x-derivative
                           Aajk = (A_next_x(j) - A_prev_x(j))/cond.STEP*2.;
                           Aakk = (A_next_x(k) - A_prev_x(k))/cond.STEP*2.;
                        }
                        // if (k == 1) // y-derivative
                        if (k == 2) { // z-derivative
                           Aajk = (A_next_z(j) - A_prev_z(j))/cond.STEP*2.;
                           Aakk = (A_next_z(k) - A_prev_z(k))/cond.STEP*2.;
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
               OrderParam<dcomplex> A_vu, A_vup, A_vum, A_vpu, A_vmu;

               // calculate the gradient at each point:
               //   loop through the OP at the point
               for (int a = 0; a < 3; a++) {    // spin index
                  for (int j = 0; j < 3; j++) { // orbital/derivative index
                     if (j == 0) temp_bc_j = Auu_BC; // choose BC for the j index
                     if (j == 1) temp_bc_j = Aww_BC;
                     if (j == 2) temp_bc_j = Avv_BC;

                     for (int k = 0; k < 3; k++) { // orbital/derivative index
                        if (k == 0) temp_bc_k = Auu_BC; // choose BC for the k index
                        if (k == 1) temp_bc_k = Aww_BC;
                        if (k == 2) temp_bc_k = Avv_BC;

                        A_vu = matrix_operator(op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size);
                        if (n_u != cond.SIZEu-1) A_vup = matrix_operator(op_vector,n_v,n_u+1,cond.SIZEv,cond.SIZEu,OP_size); // 'up' = u plus = u+1
                        if (n_u != 0) A_vum = matrix_operator(op_vector,n_v,n_u-1,cond.SIZEv,cond.SIZEu,OP_size); // 'um' = u minus = u-1
                        if (n_v != cond.SIZEv-1) A_vpu = matrix_operator(op_vector,n_v+1,n_u,cond.SIZEv,cond.SIZEu,OP_size); // 'vp' = v plus = v+1
                        if (n_v != 0) A_vmu = matrix_operator(op_vector,n_v-1,n_u,cond.SIZEv,cond.SIZEu,OP_size); // 'vm' = v minus = v-1

                        // Right now, we only use a step of h, so it is less stable...
                        //    is there a way to work around that?
                        // use BC values to calculate the gradient terms
                        if (j == 0) { // x-derivative
                           if (n_u == 0 &&            temp_bc_j.typeL == string("Neumann")) Aajj = A(a,j)/temp_bc_j.bL;
                           if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Neumann")) Aajj = A(a,j)/temp_bc_j.bR;
                           if (n_u == 0 &&            temp_bc_k.typeL == string("Neumann")) Aakj = A(a,k)/temp_bc_k.bL;
                           if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Neumann")) Aakj = A(a,k)/temp_bc_k.bR;

                           if (n_u == 0 &&            temp_bc_j.typeL == string("Dirichlet")) Aajj = ( A_vup(j) - A_vu(j) )/cond.STEP;
                           if (n_u == 0 &&            temp_bc_k.typeL == string("Dirichlet")) Aakj = ( A_vup(k) - A_vu(k) )/cond.STEP;
                           if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Dirichlet")) Aajj = ( A_vu(j) - A_vum(j) )/cond.STEP;
                           if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Dirichlet")) Aakj = ( A_vu(k) - A_vum(k) )/cond.STEP;
                        }
                        // if (j == 1) // y derivative
                        if (j == 2) { // z-derivative
                           if (n_v == 0 &&            temp_bc_j.typeB == string("Neumann")) Aajj = A(a,j)/temp_bc_j.bB;
                           if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Neumann")) Aajj = A(a,j)/temp_bc_j.bT;
                           if (n_v == 0 &&            temp_bc_k.typeB == string("Neumann")) Aakj = A(a,k)/temp_bc_k.bB;
                           if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Neumann")) Aakj = A(a,k)/temp_bc_k.bT;

                           if (n_v == 0 &&            temp_bc_j.typeB == string("Dirichlet")) Aajj = ( A_vpu(j) - A_vu(j) )/cond.STEP;
                           if (n_v == 0 &&            temp_bc_k.typeB == string("Dirichlet")) Aakj = ( A_vpu(k) - A_vu(k) )/cond.STEP;
                           if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Dirichlet")) Aajj = ( A_vu(j) - A_vmu(j) )/cond.STEP;
                           if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Dirichlet")) Aakj = ( A_vu(k) - A_vmu(k) )/cond.STEP;
                        }

                        if (k == 0) { // x-derivative
                           if (n_u == 0 &&            temp_bc_j.typeL == string("Neumann")) Aajk = A(a,j)/temp_bc_j.bL;
                           if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Neumann")) Aajk = A(a,j)/temp_bc_j.bR;
                           if (n_u == 0 &&            temp_bc_k.typeL == string("Neumann")) Aakk = A(a,k)/temp_bc_k.bL;
                           if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Neumann")) Aakk = A(a,k)/temp_bc_k.bR;

                           if (n_u == 0 &&            temp_bc_j.typeL == string("Dirichlet")) Aajk = ( A_vup(j) - A_vu(j) )/cond.STEP;
                           if (n_u == 0 &&            temp_bc_k.typeL == string("Dirichlet")) Aakk = ( A_vup(k) - A_vu(k) )/cond.STEP;
                           if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Dirichlet")) Aajk = ( A_vu(j) - A_vum(j) )/cond.STEP;
                           if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Dirichlet")) Aakk = ( A_vu(k) - A_vum(k) )/cond.STEP;
                        }
                        // if (k == 1) // y derivative
                        if (k == 2) { // z-derivative
                           if (n_v == 0 &&            temp_bc_j.typeB == string("Neumann")) Aajk = A(a,j)/temp_bc_j.bB;
                           if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Neumann")) Aajk = A(a,j)/temp_bc_j.bT;
                           if (n_v == 0 &&            temp_bc_k.typeB == string("Neumann")) Aakk = A(a,k)/temp_bc_k.bB;
                           if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Neumann")) Aakk = A(a,k)/temp_bc_k.bT;

                           if (n_v == 0 &&            temp_bc_j.typeB == string("Dirichlet")) Aajk = ( A_vpu(j) - A_vu(j) )/cond.STEP;
                           if (n_v == 0 &&            temp_bc_k.typeB == string("Dirichlet")) Aakk = ( A_vpu(k) - A_vu(k) )/cond.STEP;
                           if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Dirichlet")) Aajk = ( A_vu(j) - A_vmu(j) )/cond.STEP;
                           if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Dirichlet")) Aakk = ( A_vu(k) - A_vmu(k) )/cond.STEP;
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
            f_grad = -2./3.*(k1 + k2 + k3);
            // f_grad = gl.K1*k1 + gl.K2*k2 + gl.K3*k3; // the general way
            I_grad(ID(size,n_u,cond.SIZEu,n_v,0)) = f_grad;
            
            // this value is not normalized, while the values put into here have been...this can only be used if the GL equations we solve are not normalized
            // f_bulk = gl.B1*abs2((A*A_tran).trace()) + gl.B2*pow((A*A_dag).trace(),2) + gl.B3*(A*A_tran*A_conj*A_dag).trace() + gl.B4*(A*A_dag*A*A_dag).trace() + gl.B5*(A*A_dag*A_conj*A_tran).trace() + gl.alpha*(A*A_dag).trace();
            
            f_bulk = (gl.B1*abs2((A * A_tran).trace())
                     +gl.B2*pow( (A * A_dag).trace(), 2)
                     +gl.B3*(A * A_tran * A_conj * A_dag).trace()
                     +gl.B4*(A * A_dag * A * A_dag).trace()
                     +gl.B5*(A * A_dag * A_conj * A_tran).trace()
                     )*-1./(beta_B*9) + 2./3.*(A * A_dag).trace();
            I_bulk(ID(size,n_u,cond.SIZEu,n_v,0)) = f_bulk;

            I += ( (f_bulk + f_grad) - 1. )*cond.STEP;
         }
         I_vals(n_v) = I.real();
      }

      // save the integrand vector to plot and inspect
      WriteToFile_single_vector(I_vals,"integrand.txt");
      WriteToFile_single_vector(I_bulk,"integ_bulk.txt");
      WriteToFile_single_vector(I_grad,"integ_grad.txt");

      //cout << "The final value of f/f0 = " << integ(integ.size()-1) << endl;
      if (I.imag() >= pow(10,-8)) cout << "WARNING: imaginary part of the free-energy is not zero:" << I.imag() << endl;

      return I.real();
   }
};

#endif