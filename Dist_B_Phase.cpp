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
    Axx.bB = pow(10,-8); // a' = 0
    Axx.bT = 1.; // a = 1
    Axx.bL = pow(10,-8); // a' = 0
    Axx.bR = 1.; // a = 1
    Axx.typeB = string("Neumann");
    Axx.typeT = string("Dirichlet");
    Axx.typeL = string("Neumann");
    Axx.typeR = string("Dirichlet");

    Axz.bB = pow(10,-8); // a' = 0
    Axz.bT = 1.; // a = 1
    Axz.bL = pow(10,-8); // a' = 0
    Axz.bR = 1.; // a = 1
    Axz.typeB = string("Neumann");
    Axz.typeT = string("Dirichlet");
    Axz.typeL = string("Neumann");
    Axz.typeR = string("Dirichlet");

    Ayy.bB = pow(10,-8); // a' = 0
    Ayy.bT = 1.; // a = 1
    Ayy.bL = pow(10,-8); // a' = 0
    Ayy.bR = 1.; // a = 1
    Ayy.typeB = string("Neumann");
    Ayy.typeT = string("Dirichlet");
    Ayy.typeL = string("Neumann");
    Ayy.typeR = string("Dirichlet");

    Azx.bB = pow(10,-8); // a' = 0
    Azx.bT = 1.; // a = 1
    Azx.bL = pow(10,-8); // a' = 0
    Azx.bR = 1.; // a = 1
    Azx.typeB = string("Neumann");
    Azx.typeT = string("Dirichlet");
    Azx.typeL = string("Neumann");
    Azx.typeR = string("Dirichlet");

    Azz.bB = pow(10,-8); // a' = 0
    Azz.bT = 1.; // a = 1
    Azz.bL = 1.;  // a(L) =  1
    Azz.bR = -1.; // a(R) = -1
    Azz.typeB = string("Neumann");
    Azz.typeT = string("Dirichlet");
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
