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

    std::vector<int> no_update;
    MultiComponent_GL_Solver<VectorXcd> gls("conditions.txt");

    // start timing
    auto start = std::chrono::system_clock::now();
    gls.BuildProblem(5,Axx,Axz,Ayy,Azx,Azz);

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

    gls.Solve(guess,no_update);
    VectorXcd solution = gls.getSolution();
    gls.WriteToFile("data.txt");

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "Total time: " << elapsed.count() << " seconds." << endl;
}
