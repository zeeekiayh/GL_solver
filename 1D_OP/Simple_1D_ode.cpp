#include "includes.hpp"

//  ================
//  Global Variables
//  ================
    int SIZE;     // N; the number of points
    double ACCUR; // the minimum desired accuracy
    double STEP;  // the step size
    double B;     // the variable parameter
    int N_loop;   // count limit for iterations
    double rel_p; // relaxation parameter
//  =============

int main()
{
    std::string line;
    std::ifstream conditions ("conditions.txt");
    if (conditions.is_open())
    {
        conditions >> SIZE >> ACCUR >> STEP >> B >> N_loop >> rel_p;
        conditions.close();
    }
    else cout << "Unable to open file: conditions.txt" << endl;

    VectorXd g(SIZE); // initial guess for f
    for (int i = 0; i < SIZE; i++) g(i) = tanh(i*STEP);

    double c[6] = {(double)SIZE, ACCUR, STEP, B, (double)N_loop, rel_p};
    VectorXd f = Relaxation(c, g);

    Write_to_file(f,"points.txt");

    // free energy lost
    double Free_energy = Free_Energy(c, f);
    cout << "Free-energy result: " << Free_energy << endl;
}
