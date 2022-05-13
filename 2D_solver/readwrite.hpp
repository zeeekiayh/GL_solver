#ifndef _readwrite_hpp
#define _readwrite_hpp

#include "structures.hpp"
#include <string>
#include <sstream>

// read in the conditions from the file
inline void read_input_data(const int Nop, in_conditions & cond, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK, std::string conditions_file) {
   std::string line;
   int OP_size;
   std::ifstream conditions(conditions_file);

   // get conditions from the file
   if (conditions.is_open()) {
		conditions >> OP_size;     conditions.ignore(256,'\n');
		if(OP_size != Nop){ 
			std::cout << "wrong conditions file, OP size should be " << Nop << "; but read " << OP_size << " -- Exiting!\n"; 
			exit(1);
		}
		cond.Nop=OP_size;

		// find the line where the BOUNDARY CONDITIONS start
		while (line != "BOUNDARY CONDITIONS") {getline(conditions,line);}
		getline(conditions,line); // one extra line without data

		// Boundary conditions: >> operator for BC structure is in structures.hpp
		for(int i=0; i<Nop; i++) {
			conditions >> eta_BC[i]; conditions.ignore(256,'\n');
		}

		// find the line where the K matrix starts
		while (line != "K MATRIX") {getline(conditions,line);}
		// read the K matrix
		for (int row = 0; row < Nop; row++) {
		getline(conditions,line);

		int col=0;
		double xx,xy,yx,yy;
		char c;
		std::istringstream iss(line);
		std::string word;
		// reading iss stream word by word, words are separated by spaces
		while(iss >> word) {
			if( word[0]=='(' && word[1]==')' ) {
				gradK[row][col] << 0, 0, 0, 0; // zero matrix; 
			}
			else{
				std::istringstream word_stream(word);
				// extract data from word that has comma-separated format "(number,number,number,number)"
				word_stream >> c >> xx >> c >> xy >> c >> yx >> c >> yy >>c;
				// sscanf(word.c_str(),"(%lf,%lf,%lf,%lf)", &xx, &xy, &yx, &yy); // C-analogue: needs #include <cstdio> 
				gradK[row][col] << 	xx, xy, 
								yx, yy;
			}
			col++;
			}
		} // end for (row)

		getline(conditions,line); // empty line

		conditions >> cond.SIZEu >> cond.SIZEv; conditions.ignore(256,'\n');
		if (!cond.SIZEu) cond.SIZEu = 1; // the sizes can never be 0; set them to 1 if = 0
		if (!cond.SIZEv) cond.SIZEv = 1;

		conditions >> cond.STEP;   conditions.ignore(256,'\n');
		conditions >> cond.T;      conditions.ignore(256,'\n');
		conditions >> cond.P;      conditions.ignore(256,'\n');
		conditions >> cond.ACCUR;  conditions.ignore(256,'\n');
		conditions >> cond.N_loop; conditions.ignore(256,'\n');

		conditions.close();

		// WHY DON'T THESE ACTUALLY SET THE VALUES FOR LATER?
		// default parameters for the Convergence Accelerator
		cond.maxStore = 4;
		cond.rel_p = 0.1;
		cond.wait = 2;
   }
   else std::cout << "Unable to open file: " << conditions_file << std::endl;
}

// output the things read in so we can visually confirm the input from the conditions*.txt
inline void confirm_input_data(const int Nop, in_conditions cond, Bound_Cond eta_BC[], Eigen::Matrix2d **gradK) 
{

	for(int i=0; i<Nop; i++){
		std::cout << "bound condition OP component "<< i << " is \n" << eta_BC[i] << std::endl;
	}

	std::cout << "K-matrix:\n";
	for(int i=0; i<2*Nop; i++){ 
		for(int j=0; j<2*Nop; j++){
			std::cout << gradK[i/2][j/2](i%2,j%2)<<" "; 
			if(j%2 != 0) std::cout << "\t";
		}
		std::cout << std::endl;
		if(i%2 != 0) std::cout << "\n";
	}

	return;
}

// a simple function to map a vector's ID to it's mesh grid coordinates
// i.e., undoes: int ID(int n_u, int n_v, int n) { return n_u + Nu*n_v + grid_size*n; }
inline void unMapID(int id, int Nu, int gsize, int& n_u, int& n_v) {
	// id == (n_u + Nu*n_v + grid_size*n)
	// id % grid_size = n_u + Nu*n_v
	// (id % grid_size) % Nu = n_u
	// [(id % grid_size) - n_u] / Nu = n_v
	n_u = (id%gsize)%Nu;
	n_v = ( (id%gsize)-n_u ) / Nu;
}

// Write out the solution to a file (this is written for a 2D system)
// We can write out a single vector (like for the Free Energy) with Nop=1
inline void WriteToFile(const Eigen::VectorXcd& solution_vector, std::string file_name, const in_conditions cond) {
	std::ofstream data (file_name); // open the file for writing
	if (data.is_open()) {           // if opening was successful...

		int gsize = solution_vector.size()/cond.Nop;

		for (int i = 0; i < gsize; i++) { // loop through each grid point
			
			int Nx,Nz; // get the grid position
			unMapID(i, cond.SIZEu, gsize, Nx, Nz);

			// put the grid position first
			std::string line = std::__cxx11::to_string(cond.STEP*Nx) + std::string("\t")
							 + std::__cxx11::to_string(cond.STEP*Nz);

			for (int n = 0; n < cond.Nop; n++) { // put the element of each OP component on the same line
				line += std::string("\t") + std::__cxx11::to_string(solution_vector(i + n*gsize).real())
				      + std::string("\t") + std::__cxx11::to_string(solution_vector(i + n*gsize).imag());
			}
			data << line << std::endl; // put the whole line in the file
		}
	}
	else std::cout << "Unable to open '" << file_name << "' to write solution vector." << std::endl;
}

#endif
