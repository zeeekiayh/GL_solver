#include "complexClasses.hpp"

Three_Component_GL_Solver::Three_Component_GL_Solver(in_conditions cond) {
    this->mSize = cond.SIZEu*cond.SIZEv;
    this->OP_size = cond.OP_size;
    this->cond = cond;
}

VectorXcd Three_Component_GL_Solver::RHS() {
    Matrix3cd op, op_T, op_dag, op_conj;
    VectorXcd rhs(this->op_vector.size());

    // define an array for beta values
    double *betasNormalized;
    // ... and initialize them with the function
    betas(this->cond.T, this->cond.P, betasNormalized);

    // loop through all the OP components in the mesh
    for (int vi = 0; vi < this->OP_size; vi++) {

        for (int n_v = 0; n_v < this->cond.SIZEv; n_v++) {
        for (int n_u = 0; n_u < this->cond.SIZEu; n_u++) {
            dcomplex val; // to store the value at current mesh point

            if      (eta_BC[vi].typeB == string("Dirichlet") && n_v == 0)            val = eta_BC[vi].valB;
            else if (eta_BC[vi].typeT == string("Dirichlet") && n_v == this->cond.SIZEv-1) val = eta_BC[vi].valT;
            else
            {
                auto a = matrix_operator(this->op_vector,n_v,n_u,this->cond.SIZEv,this->cond.SIZEu,this->OP_size);
                auto axx = a(0);
                auto azz = a(2);

                // the Eta_ww_BC: perpendicular
                if (vi == 2)                 val = -1./5.*( 2.*pow(axx,2) + pow(azz,2) )*conj(azz) + 2./5.*( 2.*norm(axx) + norm(azz) )*azz + 2./5.*norm(azz)*azz - azz;
                // the Eta_uu_BC or Eta_vv_BC: parallel
                else if (vi == 0 || vi == 1) val = -1./5.*( 2.*pow(axx,2) + pow(azz,2) )*conj(axx) + 2./5.*( 2.*norm(axx) + norm(azz) )*axx + 2./5.*norm(axx)*axx - axx;

                val *= pow(this->cond.STEP,2); // because we multiplied the matrix by h^2
            }

            // insert val in rhs, in the matching location
            rhs( ID(mSize, n_u, this->cond.SIZEu, n_v, vi) ) = val;
        }
        }
    }
    return rhs;
}

double Three_Component_GL_Solver::free_energy() {
    // if (!solution.size())
    // {
    //     cout << "ERROR: cannot calculate free-energy without a solution." << endl;
    //     return 0.;
    // }

    // double beta_B = 6.*(gl.B1+gl.B2) + 2.*(gl.B3+gl.B4+gl.B5);
    // dcomplex I = 0.; // start the integral sum at 0
    // dcomplex f_bulk = 0.;
    // dcomplex f_grad = 0.;
    // dcomplex k1, k2, k3; // the 3 derivative values
    // Matrix3cd Eta_, Eta__tran, Eta__dag, Eta__conj;
    // VectorXd I_vals(cond.SIZEv); // value of the integral over distance in z
    // VectorXcd I_bulk(mSize), I_grad(mSize);
    // OrderParam<dcomplex> Eta__prev_x, Eta__next_x, Eta__prev_z, Eta__next_z; // the Eta_'s for the gradient terms
    // dcomplex Eta_ajk = 0., Eta_akj = 0., Eta_ajj = 0., Eta_akk = 0.;
    // Bound_Cond temp_bc_j, temp_bc_k;

    // // loop through all OP's in the mesh
    // for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
    //     for (int n_u = 0; n_u < cond.SIZEu; n_u++) {

    //     // calculate all the needed forms of Eta_
    //     Eta_ = matrix_operator(op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size).GetMatrixForm();
    //     Eta__tran = Eta_.transpose();
    //     Eta__conj = Eta_.conjugate();
    //     Eta__dag = Eta_.adjoint();

    //     k1 = 0., k2 = 0., k3 = 0.; // start them all at 0 for each OP
        
    //     if (n_u && n_v && n_u < cond.SIZEu-1 && n_v < cond.SIZEv-1) // if we're not at a boundary,
    //     {
    //         // get all the needed neighbors of Eta_
    //         Eta__prev_x = matrix_operator(op_vector,n_v,n_u-1,cond.SIZEv,cond.SIZEu,OP_size);
    //         Eta__next_x = matrix_operator(op_vector,n_v,n_u+1,cond.SIZEv,cond.SIZEu,OP_size);
    //         Eta__prev_z = matrix_operator(op_vector,n_v-1,n_u,cond.SIZEv,cond.SIZEu,OP_size);
    //         Eta__next_z = matrix_operator(op_vector,n_v+1,n_u,cond.SIZEv,cond.SIZEu,OP_size);

    //         // calculate the gradient at each point:
    //         //   loop through the OP at the point
    //         for (int a = 0; a < 3; a++) {       // spin index
    //             for (int j = 0; j < 3; j++) {    // orbital/derivative index
    //                 for (int k = 0; k < 3; k++) { // orbital/derivative index

    //                 // calculate the derivatives depending on the index.
    //                 // divide by h^2 later for nicer code.
    //                 if (j ==0) { // x-derivative
    //                     Eta_ajj = (Eta__next_x(j) - Eta__prev_x(j))/cond.STEP*2.;
    //                     Eta_akj = (Eta__next_x(k) - Eta__prev_x(k))/cond.STEP*2.;
    //                 }
    //                 // if (j ==1) // y-derivative
    //                 if (j ==2) { // z-derivative
    //                     Eta_ajj = (Eta__next_z(j) - Eta__prev_z(j))/cond.STEP*2.;
    //                     Eta_akj = (Eta__next_z(k) - Eta__prev_z(k))/cond.STEP*2.;
    //                 }

    //                 if (k == 0) { // x-derivative
    //                     Eta_ajk = (Eta__next_x(j) - Eta__prev_x(j))/cond.STEP*2.;
    //                     Eta_akk = (Eta__next_x(k) - Eta__prev_x(k))/cond.STEP*2.;
    //                 }
    //                 // if (k == 1) // y-derivative
    //                 if (k == 2) { // z-derivative
    //                     Eta_ajk = (Eta__next_z(j) - Eta__prev_z(j))/cond.STEP*2.;
    //                     Eta_akk = (Eta__next_z(k) - Eta__prev_z(k))/cond.STEP*2.;
    //                 }

    //                 // sum up over the indexes
    //                 k1 += norm(Eta_ajk);
    //                 k2 += conj(Eta_ajj)*Eta_akk;
    //                 k3 += conj(Eta_ajk)*Eta_akj;

    //                 // reset
    //                 Eta_ajk = 0., Eta_akj = 0., Eta_ajj = 0., Eta_akk = 0.;
    //                 } // for k
    //             } // for j
    //         } // for a
    //     } // if not @ bd
    //     else { // we're at a boundary
    //         OrderParam<dcomplex> Eta__vu, Eta__vup, Eta__vum, Eta__vpu, Eta__vmu;

    //         // calculate the gradient at each point:
    //         //   loop through the OP at the point
    //         for (int a = 0; a < 3; a++) {    // spin index
    //             for (int j = 0; j < 3; j++) { // orbital/derivative index
    //                 if (j == 0) temp_bc_j = Eta_uu_BC; // choose BC for the j index
    //                 if (j == 1) temp_bc_j = Eta_ww_BC;
    //                 if (j == 2) temp_bc_j = Eta_vv_BC;

    //                 for (int k = 0; k < 3; k++) { // orbital/derivative index
    //                 if (k == 0) temp_bc_k = Eta_uu_BC; // choose BC for the k index
    //                 if (k == 1) temp_bc_k = Eta_ww_BC;
    //                 if (k == 2) temp_bc_k = Eta_vv_BC;

    //                 Eta__vu = matrix_operator(op_vector,n_v,n_u,cond.SIZEv,cond.SIZEu,OP_size);
    //                 if (n_u != cond.SIZEu-1) Eta__vup = matrix_operator(op_vector,n_v,n_u+1,cond.SIZEv,cond.SIZEu,OP_size); // 'up' = u plus = u+1
    //                 if (n_u != 0) Eta__vum = matrix_operator(op_vector,n_v,n_u-1,cond.SIZEv,cond.SIZEu,OP_size); // 'um' = u minus = u-1
    //                 if (n_v != cond.SIZEv-1) Eta__vpu = matrix_operator(op_vector,n_v+1,n_u,cond.SIZEv,cond.SIZEu,OP_size); // 'vp' = v plus = v+1
    //                 if (n_v != 0) Eta__vmu = matrix_operator(op_vector,n_v-1,n_u,cond.SIZEv,cond.SIZEu,OP_size); // 'vm' = v minus = v-1

    //                 // Right now, we only use a step of h, so it is less stable...
    //                 //    is there a way to work around that?
    //                 // use BC values to calculate the gradient terms
    //                 if (j == 0) { // x-derivative
    //                     if (n_u == 0 &&            temp_bc_j.typeL == string("Neumann")) Eta_ajj = Eta_(a,j)/temp_bc_j.valL;
    //                     if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Neumann")) Eta_ajj = Eta_(a,j)/temp_bc_j.valR;
    //                     if (n_u == 0 &&            temp_bc_k.typeL == string("Neumann")) Eta_akj = Eta_(a,k)/temp_bc_k.valL;
    //                     if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Neumann")) Eta_akj = Eta_(a,k)/temp_bc_k.valR;

    //                     if (n_u == 0 &&            temp_bc_j.typeL == string("Dirichlet")) Eta_ajj = ( Eta__vup(j) - Eta__vu(j) )/cond.STEP;
    //                     if (n_u == 0 &&            temp_bc_k.typeL == string("Dirichlet")) Eta_akj = ( Eta__vup(k) - Eta__vu(k) )/cond.STEP;
    //                     if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Dirichlet")) Eta_ajj = ( Eta__vu(j) - Eta__vum(j) )/cond.STEP;
    //                     if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Dirichlet")) Eta_akj = ( Eta__vu(k) - Eta__vum(k) )/cond.STEP;
    //                 }
    //                 // if (j == 1) // y derivative
    //                 if (j == 2) { // z-derivative
    //                     if (n_v == 0 &&            temp_bc_j.typeB == string("Neumann")) Eta_ajj = Eta_(a,j)/temp_bc_j.valB;
    //                     if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Neumann")) Eta_ajj = Eta_(a,j)/temp_bc_j.valT;
    //                     if (n_v == 0 &&            temp_bc_k.typeB == string("Neumann")) Eta_akj = Eta_(a,k)/temp_bc_k.valB;
    //                     if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Neumann")) Eta_akj = Eta_(a,k)/temp_bc_k.valT;

    //                     if (n_v == 0 &&            temp_bc_j.typeB == string("Dirichlet")) Eta_ajj = ( Eta__vpu(j) - Eta__vu(j) )/cond.STEP;
    //                     if (n_v == 0 &&            temp_bc_k.typeB == string("Dirichlet")) Eta_akj = ( Eta__vpu(k) - Eta__vu(k) )/cond.STEP;
    //                     if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Dirichlet")) Eta_ajj = ( Eta__vu(j) - Eta__vmu(j) )/cond.STEP;
    //                     if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Dirichlet")) Eta_akj = ( Eta__vu(k) - Eta__vmu(k) )/cond.STEP;
    //                 }

    //                 if (k == 0) { // x-derivative
    //                     if (n_u == 0 &&            temp_bc_j.typeL == string("Neumann")) Eta_ajk = Eta_(a,j)/temp_bc_j.valL;
    //                     if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Neumann")) Eta_ajk = Eta_(a,j)/temp_bc_j.valR;
    //                     if (n_u == 0 &&            temp_bc_k.typeL == string("Neumann")) Eta_akk = Eta_(a,k)/temp_bc_k.valL;
    //                     if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Neumann")) Eta_akk = Eta_(a,k)/temp_bc_k.valR;

    //                     if (n_u == 0 &&            temp_bc_j.typeL == string("Dirichlet")) Eta_ajk = ( Eta__vup(j) - Eta__vu(j) )/cond.STEP;
    //                     if (n_u == 0 &&            temp_bc_k.typeL == string("Dirichlet")) Eta_akk = ( Eta__vup(k) - Eta__vu(k) )/cond.STEP;
    //                     if (n_u == cond.SIZEu-1 && temp_bc_j.typeR == string("Dirichlet")) Eta_ajk = ( Eta__vu(j) - Eta__vum(j) )/cond.STEP;
    //                     if (n_u == cond.SIZEu-1 && temp_bc_k.typeR == string("Dirichlet")) Eta_akk = ( Eta__vu(k) - Eta__vum(k) )/cond.STEP;
    //                 }
    //                 // if (k == 1) // y derivative
    //                 if (k == 2) { // z-derivative
    //                     if (n_v == 0 &&            temp_bc_j.typeB == string("Neumann")) Eta_ajk = Eta_(a,j)/temp_bc_j.valB;
    //                     if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Neumann")) Eta_ajk = Eta_(a,j)/temp_bc_j.valT;
    //                     if (n_v == 0 &&            temp_bc_k.typeB == string("Neumann")) Eta_akk = Eta_(a,k)/temp_bc_k.valB;
    //                     if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Neumann")) Eta_akk = Eta_(a,k)/temp_bc_k.valT;

    //                     if (n_v == 0 &&            temp_bc_j.typeB == string("Dirichlet")) Eta_ajk = ( Eta__vpu(j) - Eta__vu(j) )/cond.STEP;
    //                     if (n_v == 0 &&            temp_bc_k.typeB == string("Dirichlet")) Eta_akk = ( Eta__vpu(k) - Eta__vu(k) )/cond.STEP;
    //                     if (n_v == cond.SIZEv-1 && temp_bc_j.typeT == string("Dirichlet")) Eta_ajk = ( Eta__vu(j) - Eta__vmu(j) )/cond.STEP;
    //                     if (n_v == cond.SIZEv-1 && temp_bc_k.typeT == string("Dirichlet")) Eta_akk = ( Eta__vu(k) - Eta__vmu(k) )/cond.STEP;
    //                 }

    //                 k1 += norm(Eta_ajk);
    //                 k2 += conj(Eta_ajj)*Eta_akk;
    //                 k3 += conj(Eta_ajk)*Eta_akj;

    //                 // reset
    //                 Eta_ajk = 0., Eta_akj = 0., Eta_ajj = 0., Eta_akk = 0.;
    //                 } // for k
    //             } // for j
    //         } // for a
    //     }

    //     // add them all together to get the gradient term for this OP
    //     f_grad = -2./3.*(k1 + k2 + k3);
    //     // f_grad = gl.K1*k1 + gl.K2*k2 + gl.K3*k3; // the general way
    //     I_grad(ID(mSize,n_u,cond.SIZEu,n_v,0)) = f_grad;
        
    //     // this value is not normalized, while the values put into here have been...this can only be used if the GL equations we solve are not normalized
    //     // f_bulk = gl.B1*norm((Eta_*Eta__tran).trace()) + gl.B2*pow((Eta_*Eta__dag).trace(),2) + gl.B3*(Eta_*Eta__tran*Eta__conj*Eta__dag).trace() + gl.B4*(Eta_*Eta__dag*Eta_*Eta__dag).trace() + gl.B5*(Eta_*Eta__dag*Eta__conj*Eta__tran).trace() + gl.alpha*(Eta_*Eta__dag).trace();
        
    //     f_bulk = (gl.B1*norm((Eta_ * Eta__tran).trace())
    //                 +gl.B2*pow( (Eta_ * Eta__dag).trace(), 2)
    //                 +gl.B3*(Eta_ * Eta__tran * Eta__conj * Eta__dag).trace()
    //                 +gl.B4*(Eta_ * Eta__dag * Eta_ * Eta__dag).trace()
    //                 +gl.B5*(Eta_ * Eta__dag * Eta__conj * Eta__tran).trace()
    //             )*-1./(beta_B*9) + 2./3.*(Eta_ * Eta__dag).trace();
    //     I_bulk(ID(mSize,n_u,cond.SIZEu,n_v,0)) = f_bulk;

    //     I += ( (f_bulk + f_grad) - 1. )*cond.STEP;
    //     }
    //     I_vals(n_v) = I.real();
    // }

    // // save the integrand vector to plot and inspect
    // WriteToFile_single_vector(I_vals,"integrand.txt");
    // WriteToFile_single_vector(I_bulk,"integ_bulk.txt");
    // WriteToFile_single_vector(I_grad,"integ_grad.txt");

    // //cout << "The final value of f/f0 = " << integ(integ.size()-1) << endl;
    // if (I.imag() >= pow(10,-8)) cout << "WEta_RNING: imaginary part of the free-energy is not zero:" << I.imag() << endl;

    // return I.real();
    return 0.0;
}
