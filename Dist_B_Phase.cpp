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
    // these should be put into a file to read in
    Bound_Cond Axx,Axz,Ayy,Azx,Azz;
    Axx.bB = 2.;         // a' = a(B)/b
    Axx.bT = pow(10,12); // a' = 0
    Axx.bL = pow(10,12); // a' = 0
    Axx.bR = pow(10,12); // a' = 0
    Axx.typeB = string("Neumann");
    Axx.typeT = string("Neumann");
    Axx.typeL = string("Neumann");
    Axx.typeR = string("Neumann");

    Axz.bB = 2.;         // a' = a(B)/b
    Axz.bT = pow(10,12); // a' = 0
    Axz.bL = pow(10,12); // a' = 0
    Axz.bR = pow(10,12); // a' = 0
    Axz.typeB = string("Neumann");
    Axz.typeT = string("Neumann");
    Axz.typeL = string("Neumann");
    Axz.typeR = string("Neumann");

    Ayy.bB = 2.;         // a' = a(B)/b
    Ayy.bT = pow(10,12); // a' = 0
    Ayy.bL = pow(10,12); // a' = 0
    Ayy.bR = pow(10,12); // a' = 0
    Ayy.typeB = string("Neumann");
    Ayy.typeT = string("Neumann");
    Ayy.typeL = string("Neumann");
    Ayy.typeR = string("Neumann");

    Azx.bB = 2.;         // a' = a(B)/b
    Azx.bT = pow(10,12); // a' = 0
    Azx.bL = pow(10,12); // a' = 0
    Azx.bR = pow(10,12); // a' = 0
    Azx.typeB = string("Neumann");
    Azx.typeT = string("Neumann");
    Azx.typeL = string("Neumann");
    Azx.typeR = string("Neumann");

    Azz.bB = 0.;  // a = 0
    Azz.bT = pow(10,12); // a' = 0
    Azz.bL = 1.;  // a(L) =  1
    Azz.bR = -1.; // a(R) = -1
    Azz.typeB = string("Dirichlet");
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
    int comp_length = size/5; // the length of each component in guess
    int num;
    
    VectorXcd guess(size);
    for (int i = 0; i < size; i++)
    {
        num = floor(i/comp_length);
        switch (num)
        {
        case 0:
            guess(i) = 1.;
            break;
        case 1:
            guess(i) = 1.;
            break;
        case 2:
            guess(i) = 1.;
            break;
        case 3:
            guess(i) = 1.;
            break;
        case 4:
            guess(i) = ((i-4*comp_length)%c.SIZEu)/2 ? 1. : -1.; // if on the left half, put +1, otherwise put -1
            if (i < 4*comp_length + c.SIZEu) guess(i) = Azz.bB;  // set the bottom to its boundary value
            break;
        
        default:
            break;
        }
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "Problem built; time: " << elapsed.count() << " seconds." << endl;

    gls.Solve(guess,no_update);
    VectorXcd solution = gls.getSolution();
    gls.WriteToFile("data.txt");

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "Total time: " << elapsed.count() << " seconds." << endl;
}
