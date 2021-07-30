#include "includes.hpp"

//  ================
//  Global Variables
//  ================
    VectorXd guess, f;
//  ================

int main()
{
    std::vector<int> no_update;
    MultiComponent_GL_Solver<VectorXd> gls("conditions.txt","boundary_conditions.txt");
    cout << "Initialized GL Solver" << endl;

    // start timing
    auto start = std::chrono::system_clock::now();

    cout << "Building problem...";
    gls.BuildProblem();
    cout << "done" << endl;

    // make a guess
    int size = gls.getSolverMatrixSize(); // the number of columns in the matrix
    in_conditions c = gls.getConditions();
    guess.resize(size);
    auto p = gls.getConditions();
    for (int i = 0; i < size; i++) guess(i) = 1.; // i < p.SIZEu ? 0. : 1.

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

    gls.Solve(guess);

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

    VectorXd solution = gls.getSolution();
    gls.WriteToFile(solution,"data.txt");
}
