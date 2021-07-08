#include "includes.hpp"

//  ================
//  Global Variables
//  ================
    in_conditions c;
    VectorXcd guess, f;
    SparseMatrix<complex<double>> L;
    Matrix<vector<int>,3,3> var_mat;
//  ================

int main()
{
    // ReadConditions("conditions.txt", c);
    // BuildVarMatrix(var_mat);
    // // PrintVarMatrix(var_mat);
    // guess.resize(c.SIZE*2); // make some initial guess
    // for (int k = 0; k < guess.size(); k+=2)
    // {
    //     guess(k) = c.BCNp;
    //     guess(k+1) = c.BCNs;
    // }
    // BuildProblem(L, c, guess); // set up the matrix
    // f = Relaxation(L, guess, c); // use relaxation to find the solution
    // deltas f_parts = Separate(f); // separate the parallel and perpendicular
    // // save the final result to a file for Python to plot
    // Write_to_file(f_parts.para,"para_comp.txt","para_real.txt");
    // Write_to_file(f_parts.perp,"perp_comp.txt","perp_real.txt");
    // // set up for the matrix
    // vector<VectorXcd> del;
    // del.push_back(f_parts.para);
    // del.push_back(f_parts.para);
    // del.push_back(f_parts.perp);
    // vector<int> idx;
    // idx.push_back(0);
    // idx.push_back(4);
    // idx.push_back(8);
    // Matrix<VectorXcd,3,3> OP = BuildMatrix(del, idx); // the order parameter
    // // calculate the free-energy "lost"
    // complex<double> Free_energy = Free_Energy(c, OP, var_mat);
    // cout << "Free-energy result: " << Free_energy << endl;

    MultiComponent_GL_Solver<VectorXcd> gls("conditions.txt");
    gls.BuildProblem(5,Axx,Axz,Ayy,Azx,Azz);
}
