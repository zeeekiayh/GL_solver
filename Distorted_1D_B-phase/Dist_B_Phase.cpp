#include "includes.hpp"

//  ================
//  Global Variables
//  ================
    //
//  ================

int main()
{
    in_conditions c;
    std::ifstream conditions ("conditions.txt");

    if (conditions.is_open())
    {
        conditions >> c.SIZE >> c.ACCUR >> c.STEP >> c.B >> c.N_loop >> c.rel_p >> c.BC;
        conditions.close();
    }
    else cout << "Unable to open file: conditions.txt" << endl;

    // initialize all the constants
    // consts con;
    // con.At = -0.7; // actually 1/3 N_f (T/Tc - 1)
    // // B1 = N_f / (pi k_b Tc)^2 * 7 zeta(3)/240
    // con.B1 = 5; con.B2 = -con.B1/2.; con.B3 = -con.B1/2.; con.B4 = -con.B1/2.; con.B5 = -con.B1/2.;
    // // K1 = 7 zeta(3)/60 N_f xi_0^2, xi_0 = h_bar v_f / (2pi k_b Tc)
    // con.K1 = 7; con.K2 = con.K1; con.K3 = con.K1;
    // c.c = con;

    VectorXcd pa(c.SIZE), pe(c.SIZE);
    for (int k = 0; k < c.SIZE; k++)
    {
        pa(k) = 1.;
        pe(k) = 1.;
    }
    deltas del(pa, pe); // initial guess for f

    SparseMatrix<complex<double>> L_pa, L_pe;
    SetUp(L_pa, L_pe, del, c);
    deltas d = Relaxation(L_pa, L_pe, del, c);

    Write_to_file(d.para,"para_comp.txt","para_real.txt");
    Write_to_file(d.perp,"perp_comp.txt","perp_real.txt");

    // free energy "lost"
    // deltas Free_energy = Free_Energy(c, d);
    // cout << "Free-energy result: " << Free_energy << endl;
}
