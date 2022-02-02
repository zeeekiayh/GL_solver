#include "basic_gl_solver.hpp"

// =======================================
// FUNCTIONS FOR THE BASIC GL SOLVER CLASS
// CONSTRUCTORS
Basic_GL_Solver::Basic_GL_Solver() {};
Basic_GL_Solver::Basic_GL_Solver(in_conditions cond) {
	// ReadConditions(conditions_file);
	mSize = cond.SIZEu*cond.SIZEv;
	OP_size = cond.OP_size;
	this->cond = cond;
}

// MEMBER METHODS
VectorXcd Basic_GL_Solver::getSolution() const { return solution; }
int Basic_GL_Solver::getSolverMatrixSize() const { return SolverMatrix.cols(); }
SpMat_cd Basic_GL_Solver::getSolverMatrix() const { return SolverMatrix; }

// ADD CODE IN THIS FUNCTION TO BUILD ADDITIONAL D-MATRICES

// Build the derivative matrices
void Basic_GL_Solver::Build_D_Matrices() {
	// cout << "build matrices..." << endl;
	// make vectors to hold the triplets of coefficients
	vector<Tr> coeffs_u2, coeffs_w2, coeffs_uw; // coeffs_v2, coeffs_uv, coeffs_vw

	// u, v, w are orthogonal basis for the system
	//   the n's are indexes in each direction

	// Loop through the entire mesh--row/col order does not matter here
	for (int n_u = 0; n_u < cond.SIZEu; n_u++) {
		for (int n_v = 0; n_v < cond.SIZEv; n_v++) {
		int id = ID(mSize,n_u,cond.SIZEu,n_v,0);
		
		InsertCoeff_Du2(id, n_u,   n_v, -2., coeffs_u2);
		InsertCoeff_Du2(id, n_u-1, n_v,  1., coeffs_u2);
		InsertCoeff_Du2(id, n_u+1, n_v,  1., coeffs_u2);
		
		InsertCoeff_Dw2(id, n_u, n_v,  -2., coeffs_w2);
		InsertCoeff_Dw2(id, n_u, n_v-1, 1., coeffs_w2);
		InsertCoeff_Dw2(id, n_u, n_v+1, 1., coeffs_w2);
		
		InsertCoeff_Duw(id, n_u-1, n_v-1, 1./4., coeffs_uw);
		InsertCoeff_Duw(id, n_u+1, n_v-1,-1./4., coeffs_uw);
		InsertCoeff_Duw(id, n_u-1, n_v+1,-1./4., coeffs_uw);
		InsertCoeff_Duw(id, n_u+1, n_v+1, 1./4., coeffs_uw);
		}
	}

	// cout << "looped successfully" << endl;

	// initialize the D's by mSize
	Du2.resize(mSize,mSize);
	Dw2.resize(mSize,mSize);
	Duw.resize(mSize,mSize);

	// build all the D's from their coefficient triplet-vectors
	Du2.setFromTriplets(coeffs_u2.begin(), coeffs_u2.end());
	Dw2.setFromTriplets(coeffs_w2.begin(), coeffs_w2.end());
	Duw.setFromTriplets(coeffs_uw.begin(), coeffs_uw.end());

	// cout << "build matrices: done" << endl;
}

// Insert method for the Du^2 matrix derivatives
void Basic_GL_Solver::InsertCoeff_Du2(int id, int u, int v, double weight, vector<Tr>& coeffs) {
	if (u == -1 || u == cond.SIZEu);
	else coeffs.push_back(Tr(id,ID(mSize,u,cond.SIZEv,v,0),weight));
}

// insert method for the Dw^2 matrix derivatives
void Basic_GL_Solver::InsertCoeff_Dw2(int id, int u, int v, double weight, vector<Tr>& coeffs) {
	if (v == -1 || v == cond.SIZEv); // would add boundary conditions here, but
	else coeffs.push_back(Tr(id,ID(mSize,u,cond.SIZEu,v,0),weight)); // we'll use ghost points, so do nothing
}

// insert method for the mixed derivative matrices
void Basic_GL_Solver::InsertCoeff_Duw(int id, int u, int v, double weight, vector<Tr>& coeffs) {
		if (u == -1 || u == cond.SIZEu); // would add boundary conditions here,
	else if (v == -1 || v == cond.SIZEv); //  but we'll use ghost points, so do nothing
	else coeffs.push_back(Tr(id,ID(mSize,u,cond.SIZEu,v,0),weight));
}

// WRITE ADDITIONAL 'InsertCoeff_D**' METHODS HERE

// ?still needs fixed? I think it's good
// derivative matrix: 2nd-order of 1st coordinate (i.e. x)
SpMat_cd Basic_GL_Solver::Du2_BD(Bound_Cond BC, int op_elem_num) {
	// cout << "\t\t\tDu2_BD()" << endl;
	// vector<int> indexes_to_visit; // vector for debugging
	SpMat_cd Du2_copy = Du2;// the matrix that we will edit and return to not modify the original

	for (int n_v = 0; n_v < cond.SIZEv; n_v++) // loop through just the top and bottom boundary points of the mesh
	{
		// indexes for the points on the bottom side
		int id0 =         ID(mSize, 0, cond.SIZEu, n_v, 0),
			id0_connect = ID(mSize, 1, cond.SIZEu, n_v, 0);
		
		// indexes for the points on the top side
		int idN =         ID(mSize, cond.SIZEu-1, cond.SIZEu, n_v, 0),
			idN_connect = ID(mSize, cond.SIZEv-2, cond.SIZEu, n_v, 0);

		// set the values at these indexes using the ghost points,
		//   and depending on what kind of BC we have there
		if (BC.typeB == string("Neumann"))
		{
		Du2_copy.coeffRef(id0,id0) = -2. -2.*cond.STEP/BC.valB;
		Du2_copy.coeffRef(id0,id0_connect) = 2.;
		}
		else if (BC.typeB == string("Dirichlet"))
		{
		Du2_copy.coeffRef(id0,id0) = 1.;
		if (!update) no_update.push_back(ID(mSize,0,cond.SIZEu,n_v,op_elem_num));
		Du2_copy.coeffRef(id0,id0_connect) = 0.; // make sure to disconnect from the other connection
		}

		if (BC.typeT == string("Neumann"))
		{
		Du2_copy.coeffRef(idN,idN) = -2. +2.*cond.STEP/BC.valT;
		Du2_copy.coeffRef(idN,idN_connect) = 2.;
		}
		else if (BC.typeT == string("Dirichlet"))
		{
		Du2_copy.coeffRef(idN,idN) = 1.;
		if (!update) no_update.push_back(ID(mSize,cond.SIZEu-1,cond.SIZEu,n_v,op_elem_num));
		Du2_copy.coeffRef(idN,idN_connect) = 0.; // make sure to disconnect from the other connection
		}

		// // debugging (for a large matrix):
		// if (!(n_v%21-2)) indexes_to_visit.push_back(id0);
	}

	// // debugging (for a large matrix):
	// cout << endl << "ouside loop:";
	// for (auto it = indexes_to_visit.begin(); it != indexes_to_visit.end(); it++)
	// {
	//    cout << endl << "mini view:" << endl;
	//    Matrix_SubView(Du2_copy,*it-2,*it-2,7,7);
	// }

	// cout << "\t\t\tDu2_BD(): success" << endl;
	return Du2_copy;
}

// derivative matrix: 2nd-order of 3rd coordinate (i.e. z)
SpMat_cd Basic_GL_Solver::Dw2_BD(Bound_Cond BC, int op_elem_num) {
	// cout << "\t\t\tDw2_BD" << endl;
	// vector<int> indexes_to_visit; // vector for debugging
	SpMat_cd Dw2_copy = Dw2;// the matrix that we will edit and return to not modify the original

	for (int n_u = 0; n_u < cond.SIZEu; n_u++) // loop through just the top and bottom boundary points of the mesh
	{
		// indexes for the points on the bottom side
		int id0 =         ID(mSize, n_u, cond.SIZEu, 0, 0),
			id0_connect = ID(mSize, n_u, cond.SIZEu, 1, 0);
		
		// indexes for the points on the top side
		int idN =         ID(mSize, n_u, cond.SIZEu, cond.SIZEv-1, 0),
			idN_connect = ID(mSize, n_u, cond.SIZEu, cond.SIZEv-2, 0);

		// set the values at these indexes using the ghost points,
		//   and depending on what kind of BC we have there
		if (BC.typeB == string("Neumann"))
		{
		Dw2_copy.coeffRef(id0,id0) = -2. -2.*cond.STEP/BC.valB;
		Dw2_copy.coeffRef(id0,id0_connect) = 2.;
		}
		else if (BC.typeB == string("Dirichlet"))
		{
		Dw2_copy.coeffRef(id0,id0) = 1.;
		if (!update) no_update.push_back(ID(mSize,n_u,cond.SIZEu,0,op_elem_num));
		Dw2_copy.coeffRef(id0,id0_connect) = 0.; // make sure to disconnect from the other connection
		}

		if (BC.typeT == string("Neumann"))
		{
		Dw2_copy.coeffRef(idN,idN) = -2. +2.*cond.STEP/BC.valB;
		Dw2_copy.coeffRef(idN,idN_connect) = 2.;
		}
		else if (BC.typeT == string("Dirichlet"))
		{
		Dw2_copy.coeffRef(idN,idN) = 1.;
		if (!update) no_update.push_back(ID(mSize,n_u,cond.SIZEu,cond.SIZEv-1,op_elem_num));
		Dw2_copy.coeffRef(idN,idN_connect) = 0.; // make sure to disconnect from the other connection
		}

		// // debugging (for a large matrix):
		// if (!(n_u%21-2)) indexes_to_visit.push_back(id0);
	}

	// // debugging (for a large matrix):
	// cout << endl << "ouside loop:";
	// for (auto it = indexes_to_visit.begin(); it != indexes_to_visit.end(); it++)
	// {
	//    cout << endl << "mini view:" << endl;
	//    Matrix_SubView(Dw2_copy,*it-2,*it-2,7,7);
	// }

	// cout << "\t\t\tDw2_BD: success" << endl;
	return Dw2_copy;
}

// still needs fixed?
// mixed derivative: of 1st and 3rd coordinates (i.e. x & z)
SpMat_cd Basic_GL_Solver::Duw_BD(Bound_Cond BC, int op_elem_num) {
	// cout << "\t\t\tDuw_BD" << endl;
	// the matrix that we will edit and return to not modify the original
	SpMat_cd Duw_copy = Duw;

	// loop through just the boundary points of the mesh
	for (int n_v = 1; n_v < cond.SIZEv-1; n_v++) // loop through the left and right boundary points of the mesh
	{
		// indexes for the left side
		int id0 =             ID(mSize, 0, cond.SIZEu, n_v,   0), // the index of the mesh point we're at
			id0_connectB =    ID(mSize, 0, cond.SIZEu, n_v-1, 0), // index of the bottom point to connect to
			id0_connectT =    ID(mSize, 0, cond.SIZEu, n_v+1, 0), // index of the top point to connect to
			// we need to disconnect from the points that the default D matrix has
			id0_disconnectB = ID(mSize, 1, cond.SIZEu, n_v-1, 0), // index of the bottom point to disconnect from
			id0_disconnectT = ID(mSize, 1, cond.SIZEu, n_v+1, 0); // index of the top point to disconnect from
		
		// indexes for the right side
		int idN =             ID(mSize, cond.SIZEu-1, cond.SIZEu, n_v,   0), // the index of the mesh point we're at
			idN_connectB =    ID(mSize, cond.SIZEu-1, cond.SIZEu, n_v-1, 0), // index of the bottom point to connect to
			idN_connectT =    ID(mSize, cond.SIZEu-1, cond.SIZEu, n_v+1, 0), // index of the top point to connect to
			// we need to disconnect from the points that the default D matrix has
			idN_disconnectB = ID(mSize, cond.SIZEu-2, cond.SIZEu, n_v-1, 0), // index of the bottom point to disconnect from
			idN_disconnectT = ID(mSize, cond.SIZEu-2, cond.SIZEu, n_v+1, 0); // index of the top point to disconnect from
		
		// set the values at these indexes using the ghost points
		if (BC.typeL == string("Neumann"))
		{
		Duw_copy.coeffRef(id0,id0) = 0.; // disconnect from the point itself
		Duw_copy.coeffRef(id0,id0_connectT) = cond.STEP/(2.*BC.valL);
		// Duw_copy.coeffRef(id0,id0_connectT) = cond.STEP/(2.*BC.valL) (BC.valL >= pow(10,-7) ? cond.STEP/(2.*BC.valL) : 0);
		Duw_copy.coeffRef(id0,id0_connectB) = -cond.STEP/(2.*BC.valL);
		}
		else if (BC.typeL == string("Dirichlet")) Duw_copy.coeffRef(id0,id0) = 1.;
		if (BC.typeR == string("Neumann"))
		{
		Duw_copy.coeffRef(idN,idN) = 0.; // disconnect from the point itself
		Duw_copy.coeffRef(idN,idN_connectT) = cond.STEP/(2.*BC.valR);
		Duw_copy.coeffRef(idN,idN_connectB) = -cond.STEP/(2.*BC.valR);
		}
		else if (BC.typeR == string("Dirichlet")) Duw_copy.coeffRef(idN,idN) = 1.;

		// disconnect from default connections
		Duw_copy.coeffRef(id0,id0_disconnectB) = 0.;
		Duw_copy.coeffRef(id0,id0_disconnectT) = 0.;
		Duw_copy.coeffRef(idN,idN_disconnectB) = 0.;
		Duw_copy.coeffRef(idN,idN_disconnectT) = 0.;

		// If we know the VALUE at the point, add the index of it of the guess/solution
		//   vector that we already know, i.e. we don't need to update them
		if (!update && BC.typeL == string("Dirichlet")) no_update.push_back(ID(mSize,0,cond.SIZEu,n_v,op_elem_num));
		if (!update && BC.typeR == string("Dirichlet")) no_update.push_back(ID(mSize,cond.SIZEu-1,cond.SIZEu,n_v,op_elem_num));
	}

	for (int n_u = 1; n_u < cond.SIZEu-1; n_u++) // loop through the top and bottom boundary points of the mesh
	{
		// indexes for the bottom side
		int id0 =             ID(mSize, n_u,   cond.SIZEu, 0, 0), // the index of the mesh point we're at
			id0_connectL =    ID(mSize, n_u-1, cond.SIZEu, 0, 0), // index of the left point to connect to
			id0_connectR =    ID(mSize, n_u+1, cond.SIZEu, 0, 0), // index of the right point to connect to
			// we need to disconnect from the points that the default D matrix has
			id0_disconnectL = ID(mSize, n_u-1, cond.SIZEu, 1, 0), // index of the left point to disconnect from
			id0_disconnectR = ID(mSize, n_u+1, cond.SIZEu, 1, 0); // index of the right point to disconnect from
		
		// indexes for the top side
		int idN =             ID(mSize, n_u,   cond.SIZEu, cond.SIZEv-1, 0), // the index of the mesh point we're at
			idN_connectL =    ID(mSize, n_u-1, cond.SIZEu, cond.SIZEv-1, 0), // index of the left point to connect to
			idN_connectR =    ID(mSize, n_u+1, cond.SIZEu, cond.SIZEv-1, 0), // index of the right point to connect to
			// we need to disconnect from the points that the default D matrix has
			idN_disconnectL = ID(mSize, n_u-1, cond.SIZEu, cond.SIZEv-2, 0), // index of the left point to disconnect from
			idN_disconnectR = ID(mSize, n_u+1, cond.SIZEu, cond.SIZEv-2, 0); // index of the right point to disconnect from

		// set the values at these indexes using the ghost points
		if (BC.typeB == string("Neumann"))
		{
		Duw_copy.coeffRef(id0,id0) = 0.; // disconnect from the point itself
		Duw_copy.coeffRef(id0,id0_connectR) = cond.STEP/(2.*BC.valB);
		Duw_copy.coeffRef(id0,id0_connectL) = -cond.STEP/(2.*BC.valB);
		}
		else if (BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id0,id0) = 1.;
		if (BC.typeT == string("Neumann"))
		{
		Duw_copy.coeffRef(idN,idN) = 0.; // disconnect from the point itself
		Duw_copy.coeffRef(idN,idN_connectR) = cond.STEP/(2.*BC.valB);
		Duw_copy.coeffRef(idN,idN_connectL) = -cond.STEP/(2.*BC.valB);
		}
		else if (BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(idN,idN) = 1.;

		// disconnect from default connections
		Duw_copy.coeffRef(id0,id0_disconnectL) = 0.;
		Duw_copy.coeffRef(id0,id0_disconnectR) = 0.;
		Duw_copy.coeffRef(idN,idN_disconnectL) = 0.;
		Duw_copy.coeffRef(idN,idN_disconnectR) = 0.;

		// If we know the VALUE at the point, add the index of it of the guess/solution
		//   vector that we already know, i.e. we don't need to update them
		if (!update && BC.typeB == string("Dirichlet")) no_update.push_back(ID(mSize,n_u,cond.SIZEu,0,op_elem_num));
		if (!update && BC.typeT == string("Dirichlet")) no_update.push_back(ID(mSize,n_u,cond.SIZEu,cond.SIZEv-1,op_elem_num));
	}

	// special case for the corners
	int id, id_disconnect;

	// Top left
	id = ID(mSize,0,cond.SIZEu,cond.SIZEv-1,0);
	id_disconnect = ID(mSize,1,cond.SIZEu,cond.SIZEv-2,0);
	Duw_copy.coeffRef(id,id_disconnect) = 0.;
		if (BC.typeL == string("Neumann") && BC.typeT == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.valL*BC.valB);
	else if (BC.typeL == string("Dirichlet") || BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.;
		// TODO: determine if this assumption is correct, that
		//    if the function value is given for one side, we
		//    don't have to worry about the derivative condition.
		//    (x4, below)
	
	// Top right
	id = ID(mSize,cond.SIZEu-1,cond.SIZEu,cond.SIZEv-1,0);
	id_disconnect = ID(mSize,cond.SIZEu-2,cond.SIZEu,cond.SIZEv-2,0);
	Duw_copy.coeffRef(id,id_disconnect) = 0.;
		if (BC.typeR == string("Neumann") && BC.typeT == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.valR*BC.valB);
	else if (BC.typeR == string("Dirichlet") || BC.typeT == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; // here...

	// Bottom left
	id = ID(mSize,0,cond.SIZEu,0,0);
	id_disconnect = ID(mSize,1,cond.SIZEu,1,0);
	Duw_copy.coeffRef(id,id_disconnect) = 0.;
		if (BC.typeL == string("Neumann") && BC.typeB == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.valL*BC.valB);
	else if (BC.typeL == string("Dirichlet") || BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; //...here...

	// Bottom right
	id = ID(mSize,cond.SIZEu-1,cond.SIZEu,0,0);
	id_disconnect = ID(mSize,cond.SIZEu-2,cond.SIZEu,1,0);
	Duw_copy.coeffRef(id,id_disconnect) = 0.;
		if (BC.typeR == string("Neumann") && BC.typeB == string("Neumann"))     Duw_copy.coeffRef(id,id) = cond.STEP*cond.STEP/(BC.valR*BC.valB);
	else if (BC.typeR == string("Dirichlet") || BC.typeB == string("Dirichlet")) Duw_copy.coeffRef(id,id) = 1.; //...and here

	// // If we know the VALUE at the point, add the index of it of the guess/solution
	// //   vector that we already know, i.e. we don't need to update them
	// if (BC.typeL == string("Dirichlet") || BC.typeT == string("Dirichlet")) no_update.push_back(ID(sz,0,        cond.SIZEu,cond.SIZEv-1,op_elem_num));
	// if (BC.typeR == string("Dirichlet") || BC.typeT == string("Dirichlet")) no_update.push_back(ID(sz,mSize[0]-1,cond.SIZEu,cond.SIZEv-1,op_elem_num));
	// if (BC.typeL == string("Dirichlet") || BC.typeB == string("Dirichlet")) no_update.push_back(ID(sz,0,        cond.SIZEu,0,        op_elem_num));
	// if (BC.typeR == string("Dirichlet") || BC.typeB == string("Dirichlet")) no_update.push_back(ID(sz,mSize[0]-1,cond.SIZEu,0,        op_elem_num));

	// cout << "\t\t\tDuw_BD: success" << endl;
	return Duw_copy;
}

// WRITE ADDITIONAL 'D**_BD' METHODS HERE

SpMat_cd Basic_GL_Solver::BuildSolverMatrix_Test(Matrix2d **K, Bound_Cond eta_BC[]) {

	// define the solver matrix to be built
	SpMat_cd solver_mat(mSize*OP_size,mSize*OP_size);
	int x = 0, z = 1;

	// loop through the size of the OP, i.e. the size of the K matrix
	for (int m = 0; m < OP_size; m++) {
		// use the equation:
		//  [K^mn_xx D_x^2  +  K^mn_zz D_z^2  +  (K^mn_xz  +  K^mn_zx) D_xz] eta_n = f_m(eta)
		for (int n = 0; n < OP_size; n++) {
			SpMat_cd toInsert = K[m][n](x,x) * Du2_BD(eta_BC[n],m) + K[m][n](z,z) * Dw2_BD(eta_BC[n],m) + (K[m][n](x,z) + K[m][n](z,x)) * Duw_BD(eta_BC[n],m);
			solver_mat += Place_subMatrix(m, n, OP_size, toInsert);
		}
	}

	return solver_mat;
}

void Basic_GL_Solver::BuildProblem(Matrix2d **K, Bound_Cond eta_BC[]) {
	Build_D_Matrices();
	// cout << "build matrices: done" << endl;
	op_vector.resize(mSize*OP_size); // initialize OP vetor
	// cout << "op_vector resize: done" << endl;
	SolverMatrix = BuildSolverMatrix_Test(K, eta_BC);
	// SolverMatrix = BuildSolverMatrix();
	// cout << "solver matrix: done" << endl;
}

// make an educated guess using the boundary conditions
// IN PROGRESS....
VectorXcd Basic_GL_Solver::makeGuess(VectorXcd& g) {
	for (int n = 0; n < OP_size; n++) {
		if (eta_BC[n].typeB == string("Dirichlet")) {
			for (int i = cond.SIZEu*(cond.SIZEv-1); i < mSize; i++) {
				g(i) = eta_BC[n].valB;
			}
		}
		if (eta_BC[n].typeT == string("Dirichlet")) {
			for (int i = 0; i < cond.SIZEu; i++) {
				g(i) = eta_BC[n].valT;
			}
		}
		if (eta_BC[n].typeL == string("Dirichlet")) {
			for (int i = 0; i < cond.SIZEu; i++) {
				g(i) = eta_BC[n].valL;
			}
		}
		if (eta_BC[n].typeR == string("Dirichlet")) {
			for (int i = 0; i < cond.SIZEu; i++) {
				g(i) = eta_BC[n].valR;
			}
		}
	}
	
	// // eta_BC (yy) BC's
	// if (eta_BC[]vv_BC.typeB == string("Dirichlet")) {
	// 	for (int i = cond.SIZEu*(cond.SIZEv-1)+size; i < 2*size; i++) {
	// 		g(i) = eta_BC[]vv_BC.valB;
	// 	}
	// }
	// if (eta_BC[]vv_BC.typeT == string("Dirichlet")) {
	// 	for (int i = size; i < cond.SIZEu+size; i++) {
	// 		g(i) = eta_BC[]vv_BC.valT;
	// 	}
	// }
	// // eta_BC (zz) BC's
	// if (eta_BC[]ww_BC.typeB == string("Dirichlet")) {
	// 	for (int i = cond.SIZEu*(cond.SIZEv-1)+2*size; i < 3*size; i++) {
	// 		g(i) = eta_BC[]ww_BC.valB;
	// 	}
	// }
	// if (eta_BC[]ww_BC.typeT == string("Dirichlet")) {
	// 	for (int i = 2*size; i < cond.SIZEu+2*size; i++) {
	// 		g(i) = eta_BC[]ww_BC.valT;
	// 	}
	// }
	
	return g;
}