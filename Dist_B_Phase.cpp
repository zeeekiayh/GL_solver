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

    Bound_Cond Axx,Axz,Ayy,Azx,Azz;
    Axx.bB = 1.;
    Axx.bT = 1.;
    Axx.bL = 1.;
    Axx.bR = 1.;
    Axx.typeB = string("Neumann");
    Axx.typeT = string("Neumann");
    Axx.typeL = string("Neumann");
    Axx.typeR = string("Neumann");

    Axz.bB = 1.;
    Axz.bT = 1.;
    Axz.bL = 1.;
    Axz.bR = 1.;
    Axz.typeB = string("Neumann");
    Axz.typeT = string("Neumann");
    Axz.typeL = string("Neumann");
    Axz.typeR = string("Neumann");

    Ayy.bB = 1.;
    Ayy.bT = 1.;
    Ayy.bL = 1.;
    Ayy.bR = 1.;
    Ayy.typeB = string("Neumann");
    Ayy.typeT = string("Neumann");
    Ayy.typeL = string("Neumann");
    Ayy.typeR = string("Neumann");

    Azx.bB = 1.;
    Azx.bT = 1.;
    Azx.bL = 1.;
    Azx.bR = 1.;
    Azx.typeB = string("Neumann");
    Azx.typeT = string("Neumann");
    Azx.typeL = string("Neumann");
    Azx.typeR = string("Neumann");

    Azz.bB = 1.;
    Azz.bT = 1.;
    Azz.bL = 1.;
    Azz.bR = -1.;
    Azz.typeB = string("Neumann");
    Azz.typeT = string("Neumann");
    Azz.typeL = string("Dirichlet");
    Azz.typeR = string("Dirichlet");

    cout << "here" << endl;    

    MultiComponent_GL_Solver<VectorXcd> gls("conditions.txt");
    gls.BuildProblem(5,Axx,Axz,Ayy,Azx,Azz);

    cout << "here" << endl;

    // make a guess
    int size = gls.getSolverMatrixSize();
    in_conditions c = gls.getConditions();
    int midPoint = c.SIZEu/2;
    
    VectorXcd guess(size);
    for (int i = 0; i < size; i += 5)
    {
        guess(i)   = Axx.bB;
        guess(i+1) = Axx.bB;
        guess(i+2) = Ayy.bB;
        guess(i+3) = Azx.bB;
        if (i%c.SIZEu < midPoint) guess(i+4) = Azz.bL;
        else guess(i+4) = Azz.bR;
    }

    gls.Solve(guess);
    VectorXcd solution = gls.getSolution();
    // cout << "solution =\n" << solution << endl;
}
