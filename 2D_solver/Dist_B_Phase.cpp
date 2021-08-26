#include "includes.hpp"

int main()
{
    // Complex-valued GL solver
    VectorXcd guess_cplx;
    Three_Component_GL_Solver<dcomplex> gls_cplx("conditions.txt","boundary_conditions.txt");
    // cout << "Initialized complex GL Solver" << endl;

    // start timing
    auto start = std::chrono::system_clock::now();

    cout << "Building problem...";
    gls_cplx.BuildProblem();
    cout << "done" << endl;

    // make a guess
    int size = gls_cplx.getSolverMatrixSize(); // the number of columns in the matrix
    guess_cplx.resize(size);
    for (int i = 0; i < size; i++) guess_cplx(i) = 1.;

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

    gls_cplx.Solve(guess_cplx);

    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

    VectorXcd solution_cplx = gls_cplx.getSolution();
    gls_cplx.WriteToFile(solution_cplx,"cplx_data.txt");

    cout << "Free energy: " << gls_cplx.free_energy() << endl;

    return 0;
}
