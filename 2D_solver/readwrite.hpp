#ifndef _readwrite_hpp
#define _readwrite_hpp

#include "structures.hpp"
#include <string>
#include <sstream>

// read in the conditions from the file
inline void read_input_data(int& Nop, in_conditions & cond, std::vector<Bound_Cond> & eta_BC, std::string conditions_file) {
	std::string line;
	std::ifstream conditions(conditions_file);

	// get conditions from the file
	if (conditions.is_open()) {
		conditions >> Nop; conditions.ignore(256,'\n');
		cond.Nop=Nop;

		// find the line where the BOUNDARY CONDITIONS start
		while (line != "BOUNDARY CONDITIONS") {getline(conditions,line);}
		getline(conditions,line); // one extra line without data

		Bound_Cond tempBC;
		// Boundary conditions: >> operator for BC structure is in structures.hpp
		for(int i=0; i<Nop; i++) {
			conditions >> tempBC; conditions.ignore(256,'\n');
			eta_BC.push_back(tempBC);
		}

		getline(conditions,line); // empty line

		conditions >> cond.SIZEu >> cond.SIZEv; conditions.ignore(256,'\n');
		if (!cond.SIZEu) cond.SIZEu = 1; // the sizes can never be 0; set them to 1 if = 0
		if (!cond.SIZEv) cond.SIZEv = 1;

		conditions >> cond.STEP;   conditions.ignore(256,'\n');
		conditions >> cond.T;      conditions.ignore(256,'\n');
		conditions >> cond.P;      conditions.ignore(256,'\n');
		conditions >> cond.ACCUR;  conditions.ignore(256,'\n');
		conditions >> cond.N_loop; conditions.ignore(256,'\n');
		conditions >> cond.maxStore; conditions.ignore(256,'\n');
		conditions >> cond.rel_p;  conditions.ignore(256,'\n');
		conditions >> cond.wait;   conditions.ignore(256,'\n');

		conditions.close();
	}
	else std::cout << "Unable to open file: " << conditions_file << std::endl;
}

// output the things read in so we can visually confirm the input from the conditions*.txt
inline void confirm_input_data(const int Nop, in_conditions cond, std::vector<Bound_Cond> eta_BC) {
	std::cout << "Grid is " << cond.SIZEu << " x " << cond.SIZEv << std::endl;
	std::cout << "with step size h = " << cond.STEP << std::endl << std::endl;

	for(int i=0; i<Nop; i++){
		std::cout << "bound condition OP component "<< i << " is \n" << eta_BC[i] << std::endl;
	}

	return;
}

#endif