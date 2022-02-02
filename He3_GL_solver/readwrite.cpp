#include "readwrite.hpp"

// read in the conditions from the file
void read_input_data(const int Nop, in_conditions & cond, Bound_Cond eta_BC[], Matrix2d **gradK, string conditions_file) {
	string line;
	// int OP_size;
	std::ifstream conditions(conditions_file);

	// get conditions from the file
	if (conditions.is_open()) {
		conditions >> cond.OP_size; conditions.ignore(256,'\n');
		if(cond.OP_size != Nop){ 
			cout << "wrong conditions file, OP size should be " << Nop << "; but read " << cond.OP_size << " -- Exiting!\n"; 
			exit(1);
		}

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
		string word;
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
	}
	else cout << "Unable to open file: " << conditions_file << endl;
}

void confirm_input_data(const int Nop, in_conditions cond, Bound_Cond eta_BC[], Matrix2d **gradK) {

	for(int i=0; i<Nop; i++){
		cout << "bound condition OP component "<< i << " is \n" << eta_BC[i] << endl;
	}

	cout << "K-matrix:\n";
	for(int i=0; i<2*Nop; i++){ 
		for(int j=0; j<2*Nop; j++){
			cout << gradK[i/2][j/2](i%2,j%2)<<" "; 
			if(j%2 != 0) cout << "\t";
		}
		cout << endl;
		if(i%2 != 0) cout << "\n";
	}

	return;
}

void WriteToFile(VectorXcd& vec, string file_name, int OP_size, int SIZEv, int SIZEu, int STEP) {
	int mSize = SIZEv*SIZEu;

	std::ofstream data (file_name);
	if (data.is_open()) {
		for (int n_v = 0; n_v < SIZEv; n_v++) {
			for (int n_u = 0; n_u < SIZEu; n_u++) {
				string line = to_string(n_u*STEP) + string("\t")
							+ to_string(n_v*STEP) + string("\t"); // add the position components

				for (int vi = 0; vi < OP_size; vi++) {
					dcomplex element = vec(ID(mSize,n_u,SIZEu,n_v,vi));
					line += to_string(element.real())               // add the real and imaginary
						+ string("\t") + to_string(element.imag()); //  components of the 'vec' vector
					if (vi+1 < mSize) line += string("\t");
				}

				data << line << endl;
			}
		}
	}
	else cout << "Unable to open file: " << file_name << endl;
}

void WriteToFile_single_vector(VectorXd& vec, string file_name, int STEP) {
	std::ofstream data (file_name);
	if (data.is_open()) {
		for (int i = 0; i < vec.size(); i++) {
			string line = to_string(i*STEP) + string("\t") + to_string(vec(i));
			data << line << endl;
		}
	}
	else cout << "Unable to open file: " << file_name << endl;
}

void WriteToFile_single_vector(VectorXcd& vec, string file_name, int STEP) {
	std::ofstream data (file_name);
	if (data.is_open()) {
		for (int i = 0; i < vec.size(); i++) {
			string line = to_string(i*STEP) + string("\t") + to_string(vec(i).real()) + string("\t") + to_string(vec(i).imag());
			data << line << endl;
		}
	}
	else cout << "Unable to open file: " << file_name << endl;
}