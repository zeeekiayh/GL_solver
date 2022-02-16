#include "solve.hpp"

// use the relaxation method and Anderson Acceleration to solve
template <class GL_Solver>
VectorXcd Solve_Matrix_Equation(SpMat_cd& SolverMatrix, VectorXcd& guess, GL_Solver& OP, in_conditions cond, vector<int> no_update, string method) {
    
    cout << "============" << endl;
    cout << " Solving... " << endl;
    
    auto start = std::chrono::system_clock::now();

    VectorXcd f = guess, df(guess.size()); // initialize vectors

    if (no_update.size()) cout << "using no_update" << endl;

    // use LU decomposition to solve the system
    SparseLU<SparseMatrix<dcomplex>, COLAMDOrdering<int> > solver;
    solver.analyzePattern(SolverMatrix); // without this, Eigen throws: Eigen::Matrix<int, -1, 1>; ... Assertion `index >= 0 && index < size()' failed.
    solver.factorize(SolverMatrix);

    // check to see if the solver failed or not
    if (solver.info() == Eigen::Success) cout << "\tSolver: successfully built" << endl;
    else if (solver.info() == Eigen::NumericalIssue) { // for debugging non-invertable matrices
        cout << "Solver: ERROR: numerical issues... returning guess." << endl;
        // for (int k=0; k < SolverMatrix.outerSize(); ++k) for (SparseMatrix<dcomplex>::InnerIterator it(SolverMatrix,k); it; ++it) cout << "(" << it.row() << "," << it.col() << ")\t";
        return guess;
    }

    // time to prepare solving method
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "\ttime: " << elapsed.count() << " seconds." << endl << endl;

    int cts = 0; // count loops
    double err;  // to store current error
    VectorXcd rhs = OP.bulkRHS(); // the right hand side

    // the acceleration object
    converg_acceler<VectorXcd> Con_Acc(no_update); // Con_Acc(cond.maxStore,cond.wait,cond.rel_p,no_update)
    
    // loop until f converges or until it's gone too long
    do { // use relaxation

        df = solver.solve(rhs)-f; // find the change in f

        // NOTE: this "*.template next_vector<...>" may not need to be like this anymore...?
        if (method == string("acceleration")) Con_Acc.template next_vector<Matrix<dcomplex,-1,-1>>(f,df,err); // use Anderson Acceleration to converge faster
        else if (method == string("relaxation")) // use normal relaxation
        {
            f += cond.rel_p*df;
            err = df.norm()/f.norm();
        }
        else { cout << "ERROR: Unknown method type given. Returning guess." << endl; return guess; }

        // op_vector = f; // update the matrix for RHS
        rhs = OP.bulkRHS();   // update rhs
        cts++;         // increment counter

        // for debugging: to see if the solution is oscillating rather than converging
        // for (int i = 0; i < f.size(); i++) if ((i+1)%cond.SIZEu==0) cout << "\tf(" << i << ") = " << f(i) << endl;

        // output approx. percent completed                 "\033[A\33[2K\r" means clear line and restart cursor
        cout << "\033[A\33[2K\r" << "estimated: " << round((cts*100.)/cond.N_loop) << "% done" << endl;

    } while(err > cond.ACCUR && cts < cond.N_loop);

    if (err < cond.ACCUR) cout << "Found solution:" << endl;
    else cout << "Result did not converge satifactorily:" << endl;
    cout << "\titerations = " << cts << endl;
    cout << "\trelative error = " << err << endl;

    return f;
}