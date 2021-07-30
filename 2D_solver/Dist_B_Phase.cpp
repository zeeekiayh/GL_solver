#include "includes.hpp"

int main()
{
    VectorXd guess;
    std::vector<int> no_update;

    // real-valued GL solver
    Three_Component_GL_Solver<VectorXd,double> gls("conditions.txt","boundary_conditions.txt");
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

    //* Complex-valued GL solver
    VectorXcd guess_cplx;
    Three_Component_GL_Solver<VectorXcd,dcomplex> gls_cplx("conditions.txt","boundary_conditions.txt");
    cout << endl << endl << "Initialized complex GL Solver" << endl;

    // start timing
    start = std::chrono::system_clock::now();

    cout << "Building problem...";
    gls_cplx.BuildProblem();
    cout << "done" << endl;

    // make a guess
    guess_cplx.resize(size);
    p = gls_cplx.getConditions();
    for (int i = 0; i < size; i++) guess_cplx(i) = 1.; // i < p.SIZEu ? 0. : 1.

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

    gls_cplx.Solve(guess_cplx);

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

    VectorXcd solution_cplx = gls_cplx.getSolution();
    gls_cplx.WriteToFile(solution_cplx,"cplx_data.txt");
    //*/
}
