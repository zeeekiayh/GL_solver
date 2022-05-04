#include "includes.hpp"

int main()
{
    //*/ Complex-valued, 3-component GL solver
    VectorXcd guess_cplx;
    Three_Component_GL_Solver<dcomplex> gls_cplx("conditions3.txt");
    // cout << "Initialized complex GL Solver" << endl;

    // start timing
    auto start = std::chrono::system_clock::now();

    // build the problem: D matrices, solver matrix, BC's, etc.
    cout << "Building problem...";
    gls_cplx.BuildProblem();
    cout << "done" << endl;

    // make a guess
    int size = gls_cplx.getSolverMatrixSize(); // make it the size of the number of columns in the matrix
    guess_cplx.resize(size);
    for (int i = 0; i < size; i++) guess_cplx(i) = 1.; // by default, make it all 1.0

    // show the duration of setting up the problem
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

    // solve the system given the initial guess
    gls_cplx.Solve(guess_cplx);

    // display the time taken to solve
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " second(s)." << endl;

    // write the solution to file
    VectorXcd solution_cplx = gls_cplx.getSolution();
    gls_cplx.WriteToFile(solution_cplx,"cplx_data.txt");

    // calculate the free energy...need to change this
    //    to do it as a FE density over the mesh
    // cout << "Free energy: " << gls_cplx.free_energy() << endl;
    //*/

    return 0;
}
