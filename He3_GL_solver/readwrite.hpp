#ifndef _readwrite_hpp
#define _readwrite_hpp

#include "structures.hpp"
#include "includes.hpp"
#include <string>
#include <sstream>
#include <complex>
#include <eigen/Eigen/Dense>

using std::complex;
using std::ostream;
using std::string;
using namespace Eigen;

// read in the conditions from the file
void read_input_data(const int, in_conditions &, Bound_Cond*, Matrix2d**, string);

void confirm_input_data(const int, in_conditions, Bound_Cond*, Matrix2d**);

void WriteToFile(VectorXcd&, string, int, int, int, int);

void WriteToFile_single_vector(VectorXd&, string, int);

void WriteToFile_single_vector(VectorXcd&, string, int);

#endif