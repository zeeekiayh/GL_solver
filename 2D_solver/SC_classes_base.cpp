#include "SC_classes.hpp"
#include "readwrite.hpp"
#include "structures.hpp"
#include <vector>
#include <eigen/unsupported/Eigen/KroneckerProduct> // for the Kronecker product in Place_subMatrix

using namespace Eigen;
using namespace std;


// prototype for function defined in linear_eq_solver.cpp
void Solver(T_vector & f, SpMat_cd M, T_vector rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC);


typedef Triplet<double> Trpl;

// Insert the sparse matrix "sm" into the "location" (i,j) in
//   a sparse matrix of size "size" and return the matrix
SpMat_cd Place_subMatrix(int i, int j, int size, SpMat_cd sm) {
   SpMat_cd  Mij(size,size);
   Mij.insert(i,j) = 1.;
   return kroneckerProduct(Mij,sm);
}

// implementation of the derivative matrices 
// -------------------------------------------------------------------------------
// matrices for "free" space, without boundary conditions 
// -------------------------------------------------------------------------------
void SC_class :: Build_D_Matrices() {
   // make vectors to hold the triplets of coefficients
   vector<Trpl> coeffs_u2, coeffs_v2, coeffs_uv, coeffs_u, coeffs_v; // coeffs_v2, coeffs_uv, coeffs_vw
   double Wu, Wv, W; // weight for derivative coefficients 

   // u, w, v are orthogonal basis for the system 
   // we are using 2D (u, v) for now

   // Loop through the entire grid-mesh --row/col order does not matter here
   for (int u = 0; u < Nu; u++) {
      int um=(u>0) ? u-1 : u; // u-1 when exists 
      int up=(u<Nu-1) ? u+1 : u; // u+1 when exists 

      for (int v = 0; v < Nv; v++) {
         int vm=(v>0) ? v-1 : v; // v-1 when exists 
         int vp=(v<Nv-1) ? v+1 : v; // v+1 when exists 

         int id = ID(u,v,0);
         
         // when we are on the boundary, some of the indices [um,up,vm,vp] are the border points (u,v) 
         // they simply serve as indices for elements reservation of the sparse matrix that will be replaced by 
         // the appropriate boundary conditions values later 
         // Du2 -matrix triplets
         coeffs_u2.push_back( Trpl(id, ID(u ,v,0), -2.0) );
         coeffs_u2.push_back( Trpl(id, ID(um,v,0),  1.0) );
         coeffs_u2.push_back( Trpl(id, ID(up,v,0),  1.0) );
         // Du -matrix triplets
         Wu = 1.0/(up-um); // =1/2 for inside points; =1 for boundary points
         coeffs_u.push_back( Trpl(id, ID(um,v,0),  -Wu) );
         coeffs_u.push_back( Trpl(id, ID(up,v,0),   Wu) );

         // Dv2 -matrix triplets
         coeffs_v2.push_back( Trpl(id, ID(u,v ,0), -2.0) );
         coeffs_v2.push_back( Trpl(id, ID(u,vm,0),  1.0) );
         coeffs_v2.push_back( Trpl(id, ID(u,vp,0),  1.0) );
         // Dv -matrix triplets
         Wv = 1.0/(vp-vm);
         coeffs_v.push_back( Trpl(id, ID(u,vm,0),  -Wv) );
         coeffs_v.push_back( Trpl(id, ID(u,vp,0),   Wv) );

         // Duv -matrix triplets for all points on the grid
         W = Wu * Wv; // inverse area of the "grid points box" for the mixed derivative 
         // W=1/4 inside, W=1/2 on the sides, W=1 in the corner points
         coeffs_uv.push_back( Trpl(id, ID(um,vm,0),  W) );
         coeffs_uv.push_back( Trpl(id, ID(up,vp,0),  W) );
         coeffs_uv.push_back( Trpl(id, ID(um,vp,0), -W) );
         coeffs_uv.push_back( Trpl(id, ID(up,vm,0), -W) );
      }
   }


   // initialize the D's by grid_size
   Du2.resize( grid_size , grid_size );
   Dv2.resize( grid_size , grid_size );
   Duv.resize( grid_size , grid_size );
   Du.resize ( grid_size , grid_size );
   Dv.resize ( grid_size , grid_size );

   // build all the D's from their coefficient triplet-vectors
   Du2.setFromTriplets(coeffs_u2.begin(), coeffs_u2.end());
   Dv2.setFromTriplets(coeffs_v2.begin(), coeffs_v2.end());
   Duv.setFromTriplets(coeffs_uv.begin(), coeffs_uv.end());
   Du.setFromTriplets(coeffs_u.begin(), coeffs_u.end());
   Dv.setFromTriplets(coeffs_v.begin(), coeffs_v.end());
}

// -------------------------------------------------------------------------------
// add boundary conditions to the derivative matrices, by changing some of the 
// matrix elements that were pre-reserved in Build_D_Matrices() above
// -------------------------------------------------------------------------------
SpMat_cd    SC_class::Du2_BD 	 (Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC) {
   SpMat_cd Du2_copy = Du2;// the matrix that we will edit and return to not modify the original
   rhsBC = T_vector::Zero(vect_size); 

   for (int v = 0; v < Nv; v++) // go over the left and right boundaries from bottom to top 
   {
      // indexes for the points on the left side
      int id0 =         ID(0, v, 0),
          id0_connect = ID(1, v, 0);

      // indexes for the points on the right side
      int idN =         ID(Nu-1, v, 0),
          idN_connect = ID(Nu-2, v, 0);

      // set the values at these indexes using the ghost points,
      //   and adjust the RHS vector depending on what kind of BC we have there
      Du2_copy.coeffRef(id0,id0) = -2. -2.*h/BC.slipL;
      Du2_copy.coeffRef(id0,id0_connect) = 2.;
      if (BC.typeL == std::string("D"))
      {
         int id=ID(0,v,op_component);
         rhsBC(id) = -2.*h/BC.slipL*initOPvector(id);
      }

      Du2_copy.coeffRef(idN,idN) = -2. -2.*h/BC.slipR;
      Du2_copy.coeffRef(idN,idN_connect) = 2.;
      if (BC.typeR == std::string("D"))
      {
         int id=ID(Nu-1,v,op_component);
         rhsBC(id) = -2.*h/BC.slipR*initOPvector(id);
      }
   }
   return Du2_copy;
}

SpMat_cd    SC_class::Dv2_BD 	 (Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC) {
   SpMat_cd Dv2_copy = Dv2;// the matrix that we will edit and return to not modify the original
   rhsBC = T_vector::Zero(vect_size); 

   for (int u = 0; u < Nu; u++) // go over the top and bottom boundaries left to right 
   {
      // indexes for the points on the bottom side
      int id0 =         ID(u, 0, 0),
          id0_connect = ID(u, 1, 0);

      // indexes for the points on the top side
      int idN =         ID(u,Nv-1, 0),
          idN_connect = ID(u,Nv-2, 0);

      // set the values at these indexes using the ghost points,
      //   and adjust the RHS vector depending on what kind of BC we have there

      Dv2_copy.coeffRef(id0,id0) = -2. -2.*h/BC.slipB;
      Dv2_copy.coeffRef(id0,id0_connect) = 2.;
      // Dv2_copy.coeffRef(id0,id0) = (BC.typeB == std::string("N")) ? -2. -2.*h/BC.slipB : 1.;
      // Dv2_copy.coeffRef(id0,id0_connect) = (BC.typeB == std::string("N")) ? 2. : 0.;

      if (BC.typeB == std::string("D")) // Dirichlet
      {
         int id=ID(u,0,op_component);
		rhsBC(id) = -2.*h/BC.slipB*initOPvector(id); 
         	//rhsBC(id) = initOPvector(id);
      }

      Dv2_copy.coeffRef(idN,idN)= -2. -2.*h/BC.slipT;
      Dv2_copy.coeffRef(idN,idN_connect)= 2.;
      // Dv2_copy.coeffRef(idN,idN)= (BC.typeT == std::string("N")) ? -2. -2.*h/BC.slipT : 1;
      // Dv2_copy.coeffRef(idN,idN_connect)= (BC.typeT == std::string("N")) ? 2. : 0;

      if (BC.typeT == std::string("D")) // Dirichlet
      {
         int id=ID(u,Nv-1,op_component);
		rhsBC(id) = -2.*h/BC.slipT*initOPvector(id);
         //	rhsBC(id) = initOPvector(id);
      }
   }
   return Dv2_copy;
}

SpMat_cd    SC_class::Duv_BD 	 (Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC) {
   // the matrix that we will edit and return to not modify the original
   SpMat_cd Duv_copy = Duv;
   rhsBC = T_vector::Zero(vect_size); 

   // the default general formula for the mixed derivative matrix is in the Build_D_Matrices() function
   // here we only make it more accurate if we have Neumann boundary conditions 

   // we'll go around the grid boundary: RIGHT -> TOP -> LEFT -> BOTTOM 
   // we'll need to create new connections coefficients, and also eliminate some unneeded connections 
   int id, id_connectP, id_connectN, id_disconnectP, id_disconnectN, id_disconnectC; 

   // RIGHT
   if (BC.typeR == std::string("N")) // can do better than default 
   for (int v = 1; v < Nv-1; v++) 
   {
      id =             ID(Nu-1, v,   0); 
      id_connectP =    ID(Nu-1, v+1, 0); 
      id_connectN =    ID(Nu-1, v-1, 0); 
      id_disconnectP = ID(Nu-2, v+1, 0); 
      id_disconnectN = ID(Nu-2, v-1, 0); 

      Duv_copy.coeffRef(id,id_connectP) = -h/(2.*BC.slipR);
      Duv_copy.coeffRef(id,id_connectN) =  h/(2.*BC.slipR);

      // disconnect 
      Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectN) = 0.0; 
   }

   // TOP RIGHT corner 
   if (BC.typeR == std::string("N") && BC.typeT == std::string("N")){     
      id =             ID(Nu-1,Nv-1,0);
      id_disconnectP = ID(Nu-2,Nv-1,0);
      id_disconnectN = ID(Nu-1,Nv-2,0);
      id_disconnectC = ID(Nu-2,Nv-2,0);

      // if one of the slip lengths has notbeen initialized,
		//   we don't want it to blow up, so act as if it were 1
		Duv_copy.coeffRef(id,id) = h*h/(BC.slipR * BC.slipT);
      Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectN) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectC) = 0.0;
   }

   // TOP
   if (BC.typeT == std::string("N")) // can do better than default 
   for (int u = 1; u < Nu-1; u++) 
   {
      id =             ID(u,   Nv-1, 0); 
      id_connectP =    ID(u+1, Nv-1, 0); 
      id_connectN =    ID(u-1, Nv-1, 0); 
      id_disconnectP = ID(u+1, Nv-2, 0); 
      id_disconnectN = ID(u-1, Nv-2, 0); 

      Duv_copy.coeffRef(id,id_connectP) = -h/(2.*BC.slipT);
      Duv_copy.coeffRef(id,id_connectN) =  h/(2.*BC.slipT);

      // disconnect 
      Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectN) = 0.0; 
   }

   // TOP LEFT corner 
   if (BC.typeL == std::string("N") && BC.typeT == std::string("N")){     
      id =             ID(0,Nv-1,0);
      id_disconnectP = ID(1,Nv-1,0);
      id_disconnectN = ID(0,Nv-2,0);
      id_disconnectC = ID(1,Nv-2,0);

      // if one of the slip lengths has notbeen initialized,
		//   we don't want it to blow up, so act as if it were 1
		Duv_copy.coeffRef(id,id) = -h*h/(BC.slipL * BC.slipT);
      Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectN) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectC) = 0.0;
   }

   // LEFT
   if (BC.typeL == std::string("N")) // can do better than default 
   for (int v = 1; v < Nv-1; v++) 
   {
      id =             ID(0, v,   0); 
      id_connectP =    ID(0, v+1, 0); 
      id_connectN =    ID(0, v-1, 0); 
      id_disconnectP = ID(1, v+1, 0); 
      id_disconnectN = ID(1, v-1, 0); 

      Duv_copy.coeffRef(id,id_connectP) =  h/(2.*BC.slipL);
      Duv_copy.coeffRef(id,id_connectN) = -h/(2.*BC.slipL);

      // disconnect 
      Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectN) = 0.0; 
   }

   // BOTTOM LEFT corner 
   if (BC.typeL == std::string("N") && BC.typeB == std::string("N")){     
      id =             ID(0,0,0);
      id_disconnectP = ID(1,0,0);
      id_disconnectN = ID(0,1,0);
      id_disconnectC = ID(1,1,0);
      
      // if one of the slip lengths has notbeen initialized,
		//   we don't want it to blow up, so act as if it were 1
		Duv_copy.coeffRef(id,id) = h*h/(BC.slipL * BC.slipB);
      Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectN) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectC) = 0.0;
   }

   // BOTTOM
   if (BC.typeB == std::string("N")) // can do better than default 
   for (int u = 1; u < Nu-1; u++) 
   {
      id =             ID(u,   0, 0); 
      id_connectP =    ID(u+1, 0, 0); 
      id_connectN =    ID(u-1, 0, 0); 
      id_disconnectP = ID(u+1, 1, 0); 
      id_disconnectN = ID(u-1, 1, 0); 

      Duv_copy.coeffRef(id,id_connectP) =  h/(2.*BC.slipB);
      Duv_copy.coeffRef(id,id_connectN) = -h/(2.*BC.slipB);

      // disconnect 
      Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectN) = 0.0; 
   }

   // BOTTOM RIGHT corner 
   if (BC.typeR == std::string("N") && BC.typeB == std::string("N")){     
      id =             ID(Nu-1,0,0);
      id_disconnectP = ID(Nu-1,1,0);
      id_disconnectN = ID(Nu-2,0,0);
      id_disconnectC = ID(Nu-2,1,0);

      // if one of the slip lengths has notbeen initialized,
		//   we don't want it to blow up, so act as if it were 1
		Duv_copy.coeffRef(id,id) = -h*h/(BC.slipR * BC.slipB);
      Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectN) = 0.0;
      Duv_copy.coeffRef(id,id_disconnectC) = 0.0;
   }

   return Duv_copy;
}

SpMat_cd    SC_class::Du_BD 	 (Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC) {
   // the matrix that we will edit and return to not modify the original
   SpMat_cd Du_copy = Du;

   // the default general formula for the Du derivative matrix is in the Build_D_Matrices() function
   // here we only make it more accurate if we have Neumann boundary conditions 
   
   int id, id_disconnect; 

   // RIGHT
   if (BC.typeR == std::string("N")) // can do better than default 
   for (int v = 0; v < Nv; v++) 
   {
      id =             ID(Nu-1, v, 0); 
      Du_copy.coeffRef(id,id) = -1.0/(2.*BC.slipR);

      id_disconnect = ID(Nu-2, v, 0); 
      Du_copy.coeffRef(id,id_disconnect) = 0.0;
   }

   // LEFT
   if (BC.typeL == std::string("N")) // can do better than default 
   for (int v = 0; v < Nv; v++) 
   {
      id =             ID(0, v, 0); 
      Du_copy.coeffRef(id,id) = 1./(2.*BC.slipL);

      id_disconnect = ID(1, v, 0); 
      Du_copy.coeffRef(id,id_disconnect) = 0.0;
   }

   return Du_copy;
}

SpMat_cd    SC_class::Dv_BD 	 (Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC) {
   // the matrix that we will edit and return to not modify the original
   SpMat_cd Dv_copy = Dv;

   // the default general formula for the Du derivative matrix is in the Build_D_Matrices() function
   // here we only make it more accurate if we have Neumann boundary conditions 
   
   int id, id_disconnect; 

   // TOP
   if (BC.typeT == std::string("N")) // can do better than default 
   for (int u = 0; u < Nu; u++) 
   {
      id =             ID(u, Nv-1, 0); 
      Dv_copy.coeffRef(id,id) = -1.0/(2.*BC.slipT);

      id_disconnect = ID(u,Nv-2, 0); 
      Dv_copy.coeffRef(id,id_disconnect) = 0.0;
   }

   // BOTTOM
   if (BC.typeB == std::string("N")) // can do better than default 
   for (int u = 0; u < Nu; u++) 
   {
      id =             ID(u, 0, 0); 
      Dv_copy.coeffRef(id,id) = 1./(2.*BC.slipB);

      id_disconnect = ID(u, 1, 0); 
      Dv_copy.coeffRef(id,id_disconnect) = 0.0;
   }

   return Dv_copy;
}

// The method of building the solver matrix is general (at least the same for 1-, 3-, and 5-component OP's)
void SC_class :: BuildSolverMatrix( SpMat_cd & M, T_vector & rhsBC, const T_vector initOPvector, Bound_Cond eta_BC[], Matrix2d **gradK) {
   auto rhsBClocal=rhsBC;
   int x = 0, z = 1; // indexes for the K-matrix
   // Use equ. (15) in the Latex file:
   //    [K^mn_xx D_x^2  +  K^mn_zz D_z^2  +  (K^mn_xz  +  K^mn_zx) D_xz] eta_n = f_m(eta)
   // For free energy matrix we use the equation (10) and (36) in the latex file
   
   // initialize the FEgrad matrix
   FEgrad.resize(Nop*grid_size,Nop*grid_size);

   for (int m = 0; m < Nop; m++) {
      for (int n = 0; n < Nop; n++) {
         SpMat_cd toInsertM(grid_size,grid_size);
         SpMat_cd toInsertFE(grid_size,grid_size);

         if (gradK[m][n](x,x) != 0 && Nu > 1){
            toInsertM += gradK[m][n](x,x) * Du2_BD(eta_BC[n], n, initOPvector, rhsBClocal);
            if(m==n) rhsBC += gradK[m][n](x,x) * rhsBClocal;

            toInsertFE += gradK[m][n](x,x) * Du_BD(eta_BC[m], m, initOPvector, rhsBClocal).adjoint() * Du_BD(eta_BC[n], n, initOPvector, rhsBClocal);
   	   }

         if (gradK[m][n](z,z) != 0 && Nv > 1){
            toInsertM += gradK[m][n](z,z) * Dv2_BD(eta_BC[n], n, initOPvector, rhsBClocal);
            if(m==n) rhsBC += gradK[m][n](z,z) * rhsBClocal;
            
            toInsertFE += gradK[m][n](x,x) * Dv_BD(eta_BC[m], m, initOPvector, rhsBClocal).adjoint() * Dv_BD(eta_BC[n], n, initOPvector, rhsBClocal);
   	   }

         if (gradK[m][n](z,x) + gradK[m][n](x,z) != 0 && Nu > 1 && Nv > 1){
            toInsertM += (gradK[m][n](z,x) + gradK[m][n](x,z)) * Duv_BD(eta_BC[n], n, initOPvector, rhsBClocal);
            
            toInsertFE += gradK[m][n](z,x) * Dv_BD(eta_BC[m], m, initOPvector, rhsBClocal).adjoint() * Du_BD(eta_BC[n], n, initOPvector, rhsBClocal)
            		      + gradK[m][n](x,z) * Du_BD(eta_BC[m], m, initOPvector, rhsBClocal).adjoint() * Dv_BD(eta_BC[n], n, initOPvector, rhsBClocal);
   	   }

         M += Place_subMatrix( m, n, Nop, toInsertM );
         FEgrad += Place_subMatrix( m, n, Nop, toInsertFE );

	      //if(m==n) cout << "\nInserting M("<<m<<","<<n<<")= \n" << toInsertM;
      }
   }

   //cout << "M=\n" << M << endl;
   //cout << "rhsBC=\n" << rhsBC << endl;
}

void SC_class :: initialOPguess(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
	// Nop, Nu, Nv, h, ID() - are part of the SC_class, so we use them! 

	for (int n = 0; n < Nop; n++) {
		auto deltaZ = eta_BC[n].valueT - eta_BC[n].valueB;
		auto deltaX = 0.5*(eta_BC[n].valueR - eta_BC[n].valueL);
		auto middleX = 0.5*(eta_BC[n].valueR + eta_BC[n].valueL);

		// going through the entire grid
		for (int u = 0; u < Nu; u++) {
			// double x = h*u; // to put the domain wall at x=0
			double x = h*(u - Nu/2); // to put the domain wall at the middle: x=h*Nu/2
			for (int v = 0; v < Nv; v++) {
				double z = h*v;
				int id = ID(u,v,n);

				OPvector( id ) = (eta_BC[n].valueB + deltaZ * tanh(z/2)) * ( middleX + deltaX * tanh(x/2));

			}
		}
	}

   // if (debug) cout << "initial guess:\n" << OPvector << endl;

	return;
}

// a method of initializing a guess based on a previous solution
// NOTE: this funciton will mess things up if grid_size for solution is different than the grid_size for OPvector
void SC_class :: initialOPguessFromSolution(const T_vector & solution, T_vector & OPvector, std::vector<int> & no_update) {
   if (solution.size() > OPvector.size()) {
      cout << "ERROR: can't initialize a guess with a previous solution larger than the current:\n\tsolution.size() = " << solution.size() << "; OPvector.size()" << OPvector.size() << endl;
      return;
   }
   // should we fix the boundaries...? using no_update?
   for (int i = 0; i < solution.size(); i++)
      OPvector(i) = solution(i);
   for (int i = solution.size(); i < OPvector.size(); i++)
      OPvector(i) = 0.;
   // And what will happen to the solution of the new one (i.e. what OPvector will become)
   //    if the derivative matrices are different? Will it still converge?

   // "take the initial guess as 3-component solution from left to right (Dirichlet on
   //    left/right and top/bottom) - should converge within 1 iteration, since this is a solution"
}


void SC_class :: initOPguess_1DNarrowChannel(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
   // use "conditions3_1DChannel.txt"
   // cout << "initOPguess_NarrowChannel()" << endl;
   // solve for the normal 3-component system
   int Nop_init = 3;
   in_conditions cond_init;
   Bound_Cond eta_BC_init[Nop_init];
   Matrix2d **gradK_init;
   gradK_init = new Matrix2d *[Nop_init];
   for (int i = 0; i < Nop_init; i++) gradK_init[i] = new Matrix2d [Nop_init];
   read_input_data(Nop_init, cond_init, eta_BC_init, gradK_init, "conditions3_normal.txt");
   cond_init.maxStore = 5;
   cond_init.rel_p = 0.1;
   cond_init.wait = 1;
   vector<int> no_update_init;
   // for us, we need this to be the z-length of our system
   cond_init.SIZEv = Nv;
   cond_init.SIZEu = 1;
   cond_init.STEP = h;
   int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
   int VectSize_init = cond_init.Nop * GridSize_init;
   T_vector OPvector_init(VectSize_init);
   SpMat_cd M_init(VectSize_init,VectSize_init);
   T_vector rhsBC_init = T_vector::Zero(VectSize_init);
   SC_class *pSC_init = new ThreeCompHe3( Nop_init, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP );
   pSC_init->initialOPguess(eta_BC_init, OPvector_init, no_update_init);
   pSC_init->BuildSolverMatrix( M_init, rhsBC_init, OPvector_init, eta_BC_init, gradK_init );
   Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init);

   for (int n = 0; n < Nop; n++) {
      for (int v = 0; v < Nv; v++) {
         int id_forward = ID(0,v,n);
         int id_backward = ID(0,Nv-1-v,n);
         OPvector(id_forward) = OPvector_init(id_forward)*OPvector_init(id_backward);
      }
   }
	return;
}


void SC_class :: initOPguess_AzzFlip_WS2016(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
   // use "conditions5_W&S2016.txt"
   // cout << "initOPguess_AzzFlip()" << endl;
   // solve for the normal 3-component system
   int Nop_init = 3;
   in_conditions cond_init;
   Bound_Cond eta_BC_init[Nop_init];
   Matrix2d **gradK_init;
   gradK_init = new Matrix2d *[Nop_init];
   for (int i = 0; i < Nop_init; i++) gradK_init[i] = new Matrix2d [Nop_init];
   read_input_data(Nop_init, cond_init, eta_BC_init, gradK_init, "conditions3_1DChannel.txt");
   cond_init.maxStore = 5;
   cond_init.rel_p = 0.1;
   cond_init.wait = 1;
   vector<int> no_update_init;
   // for us, we need this to be the z-length of our system
   cond_init.SIZEv = Nv;
   cond_init.SIZEu = 1;
   cond_init.STEP = h;
   int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
   int VectSize_init = cond_init.Nop * GridSize_init;
   T_vector OPvector_init(VectSize_init);
   SpMat_cd M_init(VectSize_init,VectSize_init);
   T_vector rhsBC_init = T_vector::Zero(VectSize_init);
   SC_class *pSC_init = new ThreeCompHe3( Nop_init, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP );
	pSC_init->initOPguess_1DNarrowChannel(eta_BC_init, OPvector_init, no_update_init);
   pSC_init->BuildSolverMatrix( M_init, rhsBC_init, OPvector_init, eta_BC_init, gradK_init );
   Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init);

   // going through the entire grid
   for (int u = 0; u < Nu; u++) {
      for (int v = 0; v < Nv; v++) {
         for (int n = 0; n < Nop; n++) {
            double x = h*u;
				int id      = ID(u,v,n);
            int id_init = v + GridSize_init*n;

            if (n == 0 || n == 1)
               OPvector(id) = OPvector_init( id_init );
            else if (n == 2)
               OPvector(id) = tanh((x-h*Nu/2)/2) * OPvector_init( id_init );
            else if (n == 3 || n == 4)
               OPvector(id) = 0.;
			}
		}
	}
	return;
}


// "and the interesting thing to check is take the initial guess of 3-compnent solution on the
//    left side smoothly changing to 3-component solution with Azz flipped on the right side."

void SC_class :: initOPguess_AzzFlip(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
   // use "conditions5_xDeform.txt"
   // solve for the normal 3-component system for the initial guess
   int Nop_init = 3;
   in_conditions cond_init;
   Bound_Cond eta_BC_init[Nop_init];
   Matrix2d **gradK_init;
   gradK_init = new Matrix2d *[Nop_init];
   for (int i = 0; i < Nop_init; i++) gradK_init[i] = new Matrix2d [Nop_init];
   read_input_data(Nop_init, cond_init, eta_BC_init, gradK_init, "conditions3_normal.txt");
   cond_init.maxStore = 5;
   cond_init.rel_p = 0.1;
   cond_init.wait = 1;
   vector<int> no_update_init;
   // for us, we need this to be the z-length of our system
   cond_init.SIZEv = Nv;
   cond_init.SIZEu = 1;
   cond_init.STEP = h;
   int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
   int VectSize_init = cond_init.Nop * GridSize_init;
   T_vector OPvector_init(VectSize_init);
   SpMat_cd M_init(VectSize_init,VectSize_init);
   T_vector rhsBC_init = T_vector::Zero(VectSize_init);
   SC_class *pSC_init = new ThreeCompHe3( Nop_init, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP );
   pSC_init->initialOPguess(eta_BC_init, OPvector_init, no_update_init);
   pSC_init->BuildSolverMatrix( M_init, rhsBC_init, OPvector_init, eta_BC_init, gradK_init );
   Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init);

   for (int n = 0; n < Nop; n++) {
      for (int v = 0; v < Nv; v++) {
         for (int u = 0; u < Nu; u++) {
            double x = h*u;
            int id = ID(u,v,n);
            int id_init = v + GridSize_init*n;
            if (n < 3)
               OPvector(id) = OPvector_init(id_init) * ( (n==2) ? tanh((x-h*Nu/2)/2) : 1. );
            // else {
            //    auto deltaZ = eta_BC[n].valueT - eta_BC[n].valueB;
            //    auto deltaX = 0.5*(eta_BC[n].valueR - eta_BC[n].valueL);
            //    auto middleX = 0.5*(eta_BC[n].valueR + eta_BC[n].valueL);
            //    double z = h*v;
            //    OPvector( id ) = (eta_BC[n].valueB + deltaZ * tanh(z/2)) * ( middleX + deltaX * tanh(x/2));
            // }
            else OPvector(id) = 0.5;
         }
      }
   }
	return;
}


void SC_class :: initGuessWithCircularDomain(Bound_Cond eta_BC[], T_vector & OPvector, vector<int> & no_update) {
   // use "conditions5_AyyFlip.txt"

   // all boundaries should be Dirichlet => 1
   // we will create the circular domain in the middle of Ayy to be -1

   // the radius of the domain...we're assuming that the domain itself is large enough
   double r = 10;
   if (Nu*h/2 <= r || Nv*h/2 <= r)
      cout << "WARNING: the circular domain will be too large!" << endl;
   
   // location of the ceter of the grid
   int u_center = round(Nu/2);
   int v_center = round(Nv/2);

   for (int n = 0; n < Nop; n++) {
      // loop over the whole mesh
      for (int u = 0; u < Nu; u++) {
         for (int v = 0; v < Nv; v++) {
            int id = ID(u,v,n);

            if (u == 0) {
               OPvector(id) = eta_BC[n].valueL;
            } else if (u == Nu-1) {
               OPvector(id) = eta_BC[n].valueR;
            } else if (v == 0) {
               OPvector(id) = eta_BC[n].valueB;
            } else if (v == Nv-1) {
               OPvector(id) = eta_BC[n].valueT;
            } else if (n == 1) { // for only Ayy:
               OPvector(id) = tanh(  (sqrt(pow(h*(u-u_center),2) + pow(h*(v-v_center),2)) - r)/3.  );
            } else if (n == 3 || n == 4) {
               OPvector(id) = 0.;
            } else {
               OPvector(id) = 1.;
            }
         }
      }
   }
   // smooth off the initial guess a little
   for (int i = 0; i < 30; i++) {
      for (int n = 0; n < Nop; n++) {
         for (int u = 0; u < Nu; u++) {
            for (int v = 0; v < Nv; v++) {
               if ( (u > 0) && (u < Nu-1) && (v > 0) && (v < Nv-1) && (n == 0) && (n == 2) ) {
                  int id = ID(u,v,n);
                     // cout << "id = " << id << endl;
                  int idU = ID(u,v+1,n);
                     // cout << "idU = " << idU << endl;
                  int idL = ID(u-1,v,n);
                     // cout << "idL = " << idL << endl;
                  int idD = ID(u,v-1,n);
                     // cout << "idD = " << idD << endl;
                  int idR = ID(u+1,v,n);
                     // cout << "idR = " << idR << endl;

                  OPvector(id) = (OPvector(idU) + OPvector(idL) + OPvector(idD) + OPvector(idR))/4.;
                  // cout << "OPvector(id) = " << (OPvector(idU) + OPvector(idL) + OPvector(idD) + OPvector(idR))/4. << endl;
               }
            }
         }
      }
   }
}


// Some way to have this function solve for the solution to use as guess?
// void SC_class :: initialOPguessFromSolution(SC_class & SC, std::string conditions_file, Bound_Cond eta_BC[], T_vector & OPvector, std::vector<int> & no_update) {
//    // get all the information from the "conditions.txt"
// 	in_conditions cond;
// 	Bound_Cond eta_BC[Nop];      // boundary conditions for OP components
// 	Matrix2d **gradK;            // gradient coefficients in GL functional
// 	gradK = new Matrix2d *[Nop]; // the K matrix from eq. 12
// 	for (int i = 0; i < Nop; i++) gradK[i] = new Matrix2d [Nop];

// 	read_input_data(Nop, cond, eta_BC, gradK, conditions_file);
// 	//confirm_input_data(Nop, cond, eta_BC, gradK);
   
//    SC.initialOPguess(eta_BC, OPvector, no_update);
//    SC.BuildSolverMatrix()
// }

// Write out the vector to a file (this is written for a 2D system,
//    but can write 1D vectors just fine if Nv or Nu = 1).
// We can write out different kinds of vectors:
//    OPvector solution:   flag = 1; we'll use Nop from the class
//    FEvector:            flag = 0;
void SC_class :: WriteToFile(const T_vector& vector, std::string file_name, int flag) {
	std::ofstream data (file_name); // open the file for writing
	if (data.is_open()) {           // if opening was successful...
      // set precision here
      data << std::setprecision(8) << std::fixed;

      // label the columns in the output file
      data << "h*u     \th*v     ";
      if (flag == 1) { // OP vector
         // loop through all OP components...
         for (int n = 0; n < Nop; n++)
            data << "\t#" << n << ".real     #" << n << ".imag   ";
      }
      else if (flag == 0) // FE vector
         data << "\tFE";  // TODO: add more here! more columns: FE_total - FE_b, FE_total - FE_B, etc.
      data << std::endl;  // end the line

      // loop through the whole mesh...
      for (int v = 0; v < Nv; v++) {
         for (int u = 0; u < Nu; u++) {
            data << h*u << "\t" << h*v; // write the position

            if (flag == 1) { // OP vector
               // loop through all OP components...
               for (int n = 0; n < Nop; n++) {
                  int id = ID(u, v, n); // get the id
                  data << "\t" << vector(id).real() << "  " << vector(id).imag();
               }
               data << std::endl; // end the line
            } else if (flag == 0) {  // FE vector
               int id = ID(u, v, 0); // get the id
               // TODO: add more here! more columns like said above
               data << "\t" << vector(id).real() << endl; // because it should already pure real!
            }
         }
      }
	}
	else std::cout << "Unable to open '" << file_name << "' to write vector to file." << std::endl;
}
