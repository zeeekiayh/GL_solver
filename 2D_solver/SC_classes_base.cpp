#include "SC_classes.hpp"

#include "structures.hpp"
#include <vector>
#include <eigen/unsupported/Eigen/KroneckerProduct> // for the Kronecker product

using namespace Eigen;
using namespace std;

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
   // cout << "build matrices..." << endl;
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

   // cout << "looped successfully" << endl;

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

   // cout << "build matrices: done" << endl;
}

// -------------------------------------------------------------------------------
// add boundary conditions to the derivative matrices, by changing some of the 
// matrix elements that were pre-reserved in Build_D_Matrices() above
// -------------------------------------------------------------------------------
SpMat_cd    SC_class::Du2_BD 	 (Bound_Cond BC, int op_component, /*const*/ VectorXcd initOPvector, VectorXcd & rhsBC, double K_tilde)
{
   SpMat_cd Du2_copy = Du2;// the matrix that we will edit and return to not modify the original

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
			// initOPvector(id) = 0.;
         rhsBC(id) = -2.*h/BC.valueL*initOPvector(id) * K_tilde;
      }

      Du2_copy.coeffRef(idN,idN) = -2. -2.*h/BC.slipR;
      Du2_copy.coeffRef(idN,idN_connect) = 2.;
      if (BC.typeR == std::string("D"))
      {
         int id=ID(Nu-1,v,op_component);
			// initOPvector(id) = 0.;
         rhsBC(id) = -2.*h/BC.slipR*initOPvector(id) * K_tilde; // K_tilde from equation 27
      }
   }
   return Du2_copy;
}

SpMat_cd    SC_class::Dv2_BD 	 (Bound_Cond BC, int op_component, /*const*/ VectorXcd initOPvector, VectorXcd & rhsBC, double K_tilde)
{
   SpMat_cd Dv2_copy = Dv2;// the matrix that we will edit and return to not modify the original

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
      // Dv2_copy.coeffRef(id0,id0) = (BC.typeB == std::string("N")) ? -2. -2.*h/BC.slipB : 1.;

      Dv2_copy.coeffRef(id0,id0_connect) = 2.;
      // Dv2_copy.coeffRef(id0,id0_connect) = (BC.typeB == std::string("N")) ? 2. : 0.;

      if (BC.typeB == std::string("D")) // Dirichlet
      {
         int id=ID(u,0,op_component);
			// initOPvector(id) = 0.;
         // rhsBC(id) = BC.valueB;
			rhsBC(id) = -2.*h/BC.slipB*initOPvector(id) * K_tilde; // K_tilde from equation 27
      }

      Dv2_copy.coeffRef(idN,idN)= -2. -2.*h/BC.slipT;
      // Dv2_copy.coeffRef(idN,idN)= (BC.typeT == std::string("N")) ? -2. -2.*h/BC.slipT : 1;

      Dv2_copy.coeffRef(idN,idN_connect)= 2.;
      // Dv2_copy.coeffRef(idN,idN_connect)= (BC.typeT == std::string("N")) ? 2. : 0;

      if (BC.typeT == std::string("D")) // Dirichlet
      {
         int id=ID(u,Nv-1,op_component);
			// initOPvector(id) = 0.;
         // rhsBC(id) = BC.valueT;
			rhsBC(id) = -2.*h/BC.slipT*initOPvector(id) * K_tilde;
      }
   }
   return Dv2_copy;
}

SpMat_cd    SC_class::Duv_BD 	 (Bound_Cond BC, int op_component, /*const*/ VectorXcd initOPvector, VectorXcd & rhsBC)
{
   // the matrix that we will edit and return to not modify the original
   SpMat_cd Duv_copy = Duv;

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

SpMat_cd    SC_class::Du_BD 	 (Bound_Cond BC, int op_component, const VectorXcd initOPvector, VectorXcd & rhsBC)
{
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
      Du_copy.coeffRef(id,id) = -1.0/(2.*BC.valueR);

      id_disconnect = ID(Nu-2, v, 0); 
      Du_copy.coeffRef(id,id_disconnect) = 0.0;
   }

   // LEFT
   if (BC.typeL == std::string("N")) // can do better than default 
   for (int v = 0; v < Nv; v++) 
   {
      id =             ID(0, v, 0); 
      Du_copy.coeffRef(id,id) = 1./(2.*BC.valueL);

      id_disconnect = ID(1, v, 0); 
      Du_copy.coeffRef(id,id_disconnect) = 0.0;
   }

   return Du_copy;
}

SpMat_cd    SC_class::Dv_BD 	 (Bound_Cond BC, int op_component, const VectorXcd initOPvector, VectorXcd & rhsBC)
{
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
      Dv_copy.coeffRef(id,id) = -1.0/(2.*BC.valueT);

      id_disconnect = ID(u,Nv-2, 0); 
      Dv_copy.coeffRef(id,id_disconnect) = 0.0;
   }

   // BOTTOM
   if (BC.typeB == std::string("N")) // can do better than default 
   for (int u = 0; u < Nu; u++) 
   {
      id =             ID(u, 0, 0); 
      Dv_copy.coeffRef(id,id) = 1./(2.*BC.valueB);

      id_disconnect = ID(u, 1, 0); 
      Dv_copy.coeffRef(id,id_disconnect) = 0.0;
   }

   return Dv_copy;
}

// The method of building the solver matrix is general (at least the same for 1-, 3-, and 5-component OP's)
void SC_class :: BuildSolverMatrix( SpMat_cd & M, VectorXcd & rhsBC, const VectorXcd initOPvector, Bound_Cond eta_BC[], Matrix2d **gradK) {
   int x = 0, z = 1; // indexes for the K-matrix
   // Use equ. (15) in the Latex file:
   //    [K^mn_xx D_x^2  +  K^mn_zz D_z^2  +  (K^mn_xz  +  K^mn_zx) D_xz] eta_n = f_m(eta)
   for (int m = 0; m < Nop; m++) {
      for (int n = 0; n < Nop; n++) {
         SpMat_cd toInsert(grid_size,grid_size);

         if (gradK[m][n](x,x) != 0 && Nu > 1)
            toInsert += gradK[m][n](x,x) * Du2_BD(eta_BC[n], m, initOPvector, rhsBC);

         if (gradK[m][n](z,z) != 0 && Nv > 1)
            toInsert += gradK[m][n](z,z) * Dv2_BD(eta_BC[n], m, initOPvector, rhsBC, gradK[m][n](z,z));

         if (gradK[m][n](z,x) + gradK[m][n](x,z) != 0 && Nu > 1 && Nv > 1)
            toInsert += (gradK[m][n](z,x) + gradK[m][n](x,z)) * Du2_BD(eta_BC[n], m, initOPvector, rhsBC);

         M += Place_subMatrix( m, n, Nop, toInsert );
      }
   }
}

void SC_class :: initialOPguess(Bound_Cond eta_BC[], VectorXcd & OPvector, vector<int> & no_update) {
	// Nop, Nu, Nv, h, ID() - are part of the SC_class, so we use them! 

	// cout << "In 'initializeOPguess'" << endl;
	for (int n = 0; n < Nop; n++) {
		// cout << "\tn = " << n << endl;
		
		complex<double> deltaZ = eta_BC[n].valueT - eta_BC[n].valueB;
		complex<double> deltaX = eta_BC[n].valueR - eta_BC[n].valueL;

		// going through the entire grid
		for (int u = 0; u < Nu; u++) {
			double x = h*(u - Nu/2);
			for (int v = 0; v < Nv; v++) {
				double z = h*v;
				int id = ID(u,v,n);

				OPvector( id ) = ( eta_BC[n].valueB + deltaZ * tanh(z/2.0) )
									 * ( eta_BC[n].valueL + deltaX * tanh(x/2.0) );
                            // = 1.;

			}
		}
	}

	// cout << "initial guess:\n" << OPvector << endl;

	return;
}