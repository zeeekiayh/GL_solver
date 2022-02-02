#include "includes.hpp"

// returns the unique id corresponding to each op-component in the mesh.
//    n_u: mesh index along u
//    n_v: mesh index along v
//    i:   vector (OP) index
//    gsize: grid size = num elements in mesh
int ID(int gsize, int n_u, int n_u_max, int n_v, int i) { return gsize*i + n_u_max*n_v + n_u; }


// Insert the matrix "spMat" into the "location" (i,j) in
//   a sparse matrix of size "size" and return the matrix
SpMat_cd Place_subMatrix(int i, int j, int size, SpMat_cd spMat) {
	SpMat_cd Mij(size,size);
	Mij.insert(i,j) = 1.;
	return kroneckerProduct(Mij,spMat);
}


// A special function to access OP elements from a vector that
//   has been flattened from a matrix of size 'sizeU' X 'sizeV'
VectorXcd matrix_operator(VectorXcd& vec, int v, int u, int sizeV, int sizeU, int OPsize, bool debug) {
	VectorXcd op(OPsize); // initialize the OP to return
	if (debug) cout << "size = " << sizeV*sizeU << endl;
	for (int vi = 0; vi < OPsize; vi++) {
		if (debug) cout << "ID = " << ID(sizeU*sizeV,u,sizeU,v,vi) << endl;
		op(vi) = vec(ID(sizeU*sizeV,u,sizeU,v,vi)); // insert all OP components into the OrderParam object
	}
	return op;
}

// Now include the real and complexs class header files,
//   since the parent class are defined
// #include "complexClasses.hpp"