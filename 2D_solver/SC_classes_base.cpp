#include "SC_classes.hpp"
#include "readwrite.hpp"
#include "structures.hpp"
#include <eigen/unsupported/Eigen/KroneckerProduct> // for the Kronecker product in Place_subMatrix

using namespace Eigen;
using namespace std;

// prototype for function defined in linear_eq_solver.cpp
void Solver(T_vector & f, const SpMat_cd & M, T_vector rhsBC, in_conditions cond, vector<int> no_update, SC_class *SC);

typedef Triplet<double> Trpl;

// Insert the sparse matrix "sm" into the "location" (i,j) in a sparse matrix of size "size" and return the matrix
SpMat_cd Place_subMatrix(int i, int j, int size, SpMat_cd sm) {
   SpMat_cd  Mij(size,size);
   Mij.insert(i,j) = 1.;
   return kroneckerProduct(Mij,sm);
}

// implementation of the derivative matrices 
// -------------------------------------------------------------------------------
   // matrices for "free" space: without boundary conditions
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
            Wu = (up-um) ? 1.0/(up-um) : 0.0; // =1/2 for inside points; =1 for boundary points, 0 if no range
            coeffs_u.push_back( Trpl(id, ID(um,v,0),  -Wu) );
            coeffs_u.push_back( Trpl(id, ID(up,v,0),   Wu) );

            // Dv2 -matrix triplets
            coeffs_v2.push_back( Trpl(id, ID(u,v ,0), -2.0) );
            coeffs_v2.push_back( Trpl(id, ID(u,vm,0),  1.0) );
            coeffs_v2.push_back( Trpl(id, ID(u,vp,0),  1.0) );
            // Dv -matrix triplets
            Wv = (vp-vm) ? 1.0/(vp-vm) : 0.0; 
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

      return;
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

         Du2_copy.coeffRef(idN,idN) = -2. -2.*h/BC.slipR;
         Du2_copy.coeffRef(idN,idN_connect) = 2.;
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

         Dv2_copy.coeffRef(idN,idN)= -2. -2.*h/BC.slipT;
         Dv2_copy.coeffRef(idN,idN_connect)= 2.;
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

      Duv_copy.coeffRef(id,id) = -h*h/(BC.slipR * BC.slipB);
         Duv_copy.coeffRef(id,id_disconnectP) = 0.0;
         Duv_copy.coeffRef(id,id_disconnectN) = 0.0;
         Duv_copy.coeffRef(id,id_disconnectC) = 0.0;
      }

      return Duv_copy;
   }

   SpMat_cd    SC_class::Du_BD 	 (Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC) {
      SpMat_cd Du_copy = Du; // the matrix that we will edit and return to not modify the original
      // the default general formula for the Du derivative matrix is in the Build_D_Matrices() function
      // here we only make it more accurate if we have Neumann boundary conditions 
      
      int id, id_disconnect;
      // RIGHT
      if (BC.typeR == std::string("N")) // can do better than default 
      for (int v = 0; v < Nv; v++) 
      {
         id =             ID(Nu-1, v, 0); 
         Du_copy.coeffRef(id,id) = -h/(BC.slipR);

         id_disconnect = ID(Nu-2, v, 0); 
         Du_copy.coeffRef(id,id_disconnect) = 0.0;
      }

      // LEFT
      if (BC.typeL == std::string("N")) // can do better than default 
      for (int v = 0; v < Nv; v++) 
      {
         id =             ID(0, v, 0); 
         Du_copy.coeffRef(id,id) = h/(BC.slipL);

         id_disconnect = ID(1, v, 0); 
         Du_copy.coeffRef(id,id_disconnect) = 0.0;
      }

      return Du_copy;
   }

   SpMat_cd    SC_class::Dv_BD 	 (Bound_Cond BC, int op_component, const T_vector & initOPvector, T_vector & rhsBC) {
      SpMat_cd Dv_copy = Dv; // the matrix that we will edit and return to not modify the original
      // the default general formula for the Du derivative matrix is in the Build_D_Matrices() function
      // here we only make it more accurate if we have Neumann boundary conditions 
      
      int id, id_disconnect;
      // TOP
      if (BC.typeT == std::string("N")) // can do better than default 
      for (int u = 0; u < Nu; u++) 
      {
         id =             ID(u, Nv-1, 0); 
         Dv_copy.coeffRef(id,id) = -h/(BC.slipT);

         id_disconnect = ID(u,Nv-2, 0); 
         Dv_copy.coeffRef(id,id_disconnect) = 0.0;
      }

      // BOTTOM
      if (BC.typeB == std::string("N")) // can do better than default 
      for (int u = 0; u < Nu; u++) 
      {
         id =             ID(u, 0, 0); 
         Dv_copy.coeffRef(id,id) = h/(BC.slipB);

         id_disconnect = ID(u, 1, 0); 
         Dv_copy.coeffRef(id,id_disconnect) = 0.0;
      }

      return Dv_copy;
   }
// end of derivative matrices
// -------------------------------------------------------------------------------

T_scalar profile_dbl_tanh(double x, double x0, T_scalar v1, T_scalar v2) { return (v1+(1.0-v1)*tanh(x/2.0)) * (v2+(1.0-v2)*tanh((x0-x)/2.0)); };
T_scalar profile_dom_wall(double x, double x0, T_scalar v1, T_scalar v2) { return 0.5*(v2+v1) + 0.5*(v2-v1)*tanh((x-x0)/2.0); };

void SC_class :: initialOPguess(std::vector<Bound_Cond> eta_BC, T_vector & OPvector, vector<int> & no_update) {
	// Nop, Nu, Nv, h, ID() - are part of the SC_class, so we use them! 
	T_scalar (* profileX)(double, double, T_scalar, T_scalar); 
	T_scalar (* profileZ)(double, double, T_scalar, T_scalar); 
	double x0, z0;

	for (int n = 0; n < Nop; n++) {
		// from boundary conditions determine which profile function to use 
		if ( real( eta_BC[n].valueL * eta_BC[n].valueR) < -0.1 && Nu>1) {
			x0 = (Nu/2)*h; 
			profileX = profile_dom_wall; 
		}else{
			x0 = (Nu-1)*h; 
			profileX = profile_dbl_tanh; 
		}

		if ( real( eta_BC[n].valueT * eta_BC[n].valueB) < -0.1 && Nv>1 ) {
			z0 = (Nv/2)*h; 
			profileZ = profile_dom_wall; 
		}else{
			z0 = (Nv-1)*h; 
			profileZ = profile_dbl_tanh; 
		}
		// we set overall amplitude to zero if OP[n]=0 on all boundaries  
		double C = abs(eta_BC[n].valueL)+abs(eta_BC[n].valueB)+abs(eta_BC[n].valueR)+abs(eta_BC[n].valueT); 
		if (C>0.01) C=1; else C=0;

		// going through the entire grid
		for (int u = 0; u < Nu; u++) { double x = h*u; 
			for (int v = 0; v < Nv; v++) { double z = h*v;  
				int id = ID(u,v,n);

				OPvector( id ) = C * profileX(x, x0, eta_BC[n].valueL, eta_BC[n].valueR ) 
						   * profileZ(z, z0, eta_BC[n].valueB, eta_BC[n].valueT );

				if(u==0  && eta_BC[n].typeL=="D" && Nu>1) no_update.push_back( id );
				if(u==Nu-1 && eta_BC[n].typeR=="D" && Nu>1) no_update.push_back( id );
				if(v==0  && eta_BC[n].typeB=="D" && Nv>1) no_update.push_back( id );
				if(v==Nv-1 && eta_BC[n].typeT=="D" && Nv>1) no_update.push_back( id );
			}
		}
	}

	return;
}

// a method of initializing a guess based on a previous solution
// Works, with minimal changes for slab or semi-infinite system. This is determined by the boundary conditions 
void SC_class :: initOPguess_special(in_conditions cond, std::vector<Bound_Cond> eta_BC, T_vector & OPvector, Eigen::VectorXd & FE_ref, std::vector<int> & no_update) {
   vector<int> no_update_init;
   in_conditions cond_init=cond; // copy the conditions, 
   // but change them appropriately to fit this guess 
      cond_init.Nop=3;
      cond_init.SIZEu = 1; // make it vertical, effectively setting u=0 
   vector<Bound_Cond> eta_BC_init(Nop); // change the Top boundary conditons 
   	eta_BC_init[0]=eta_BC[0];
   	eta_BC_init[1]=eta_BC[1];
   	eta_BC_init[2]=eta_BC[2];
	// if Azz > 0 on top we want semi-infinite system along z, and we improve convergence of initial solution by setting D on top
	if(abs(eta_BC[2].valueT)>0.1) { 
		eta_BC_init[0].typeT="D"; eta_BC_init[0].slipT=1e-10; 
		eta_BC_init[1].typeT="D"; eta_BC_init[1].slipT=1e-10;
		eta_BC_init[2].typeT="D"; eta_BC_init[2].slipT=1e-10;
	}
   int GridSize_init = cond_init.SIZEu * cond_init.SIZEv;
   int VectSize_init = cond_init.Nop * GridSize_init;
   T_vector OPvector_init(VectSize_init);
   SpMat_cd M_init(VectSize_init,VectSize_init);
   T_vector rhsBC_init = T_vector::Zero(VectSize_init);
   VectorXd freeEb_init(GridSize_init), freeEg_init(GridSize_init), FEdens_init(GridSize_init); // free energy on the grid points
   VectorXd FEref_init=VectorXd::Constant(GridSize_init,-1.0);  // reference free energy, -1=bulk uniform value by default
   FEdens_init = FEref_init; // FE density on input is the uniform bulk free energy

   SC_class *pSC_init = new ThreeCompHe3( cond_init.Nop, cond_init.SIZEu, cond_init.SIZEv, cond_init.STEP );
   pSC_init->initialOPguess(eta_BC_init, OPvector_init, no_update_init);
   pSC_init->BuildSolverMatrix(M_init, rhsBC_init, OPvector_init, eta_BC_init);

   Solver(OPvector_init, M_init, rhsBC_init, cond_init, no_update_init, pSC_init);
   double totalFE = pSC_init->FreeEn(OPvector_init, cond_init, FEdens_init, freeEb_init, freeEg_init); 
   cout << "initial 1D Guess total FE energy relative to bulk's is " << totalFE << "\n";
   pSC_init->WriteAllToFile(OPvector_init, FEdens_init, freeEb_init, FEref_init, "initial_guess_OP"+to_string(cond_init.Nop)+".txt");

	double R=15.0;

	for (int v = 0; v < Nv; v++) {
		double z=h*v;
	   for (int u = 0; u < Nu; u++) {
		double x = h*(u-Nu/2);
		double rho=sqrt(x*x+z*z);

		for (int n = 0; n < Nop; n++) {
			int id = ID(u,v,n);
			int id_init = v + GridSize_init*n;

			if     (n==0) OPvector(id) = OPvector_init(id_init);  
			else if(n==1) OPvector(id) = OPvector_init(id_init);  
			else if(n==2) OPvector(id) = OPvector_init(id_init) * tanh((R-rho)/2); // semi-circle domain in Azz near a wall
			//else if(n==2) OPvector(id) = OPvector_init(id_init) * tanh(-x/2); // domain wall in Azz component along x
			else if(n==3) OPvector(id) = 0.0/cosh(x/8)/cosh(z/8);
		   	else 		  OPvector(id) = 0.0;

			// Don't need to update values that we already know from the BC's
			if(u==0  && eta_BC[n].typeL=="D" && Nv>1) no_update.push_back( id );
			if(u==Nu-1 && eta_BC[n].typeR=="D" && Nv>1) no_update.push_back( id );
			if(v==0  && eta_BC[n].typeB=="D" && Nu>1) no_update.push_back( id );
			if(v==Nv-1 && eta_BC[n].typeT=="D" && Nu>1) no_update.push_back( id );
			if( n==2 &&  abs(rho-R) < 0.4*h) no_update.push_back( id ); // to fix the radius of the semi-circle domain
		}
		FE_ref( ID(u,v,0) ) = FEdens_init( v );

	   }
	}
	delete pSC_init;
	return;
}

void SC_class :: WriteAllToFile(const T_vector& solution, const Eigen::VectorXd & FEdens, const Eigen::VectorXd & FEbulk, const Eigen::VectorXd & FEref, std::string file_name)
{
	std::ofstream data (file_name); // open the file for writing
	if (data.is_open()) {           // if opening was successful...
      // set precision here
      data << std::setprecision(8) << std::fixed;

      // label the columns in the output file
      data << "#x/xi     \t\t z/xi     ";

      // loop through all OP components...
      for (int n = 0; n < Nop; n++)
         data << "\tRe(OP" << n << ") \tIm(OP" << n << ")  ";

      // add more here if needed 
      data << "\tFEdensity  \tbulk_FE   \tgrad_FE \t FE-FEref";
      data << std::endl; 

	std::ofstream fexint ("free_energy_x_integrated.dat"); 
      fexint << std::setprecision(8) << std::fixed << "#z/xi   \t FE_xintegrated\n";
      // loop through the whole mesh...
      double x_shift=(Nu/2)*h;
      double z_shift=0.0;//(Nv/2)*h; 
      for (int v = 0; v < Nv; v++) {
		double FE_xintegrated=0.0;
         for (int u = 0; u < Nu; u++) {
            data << h*u-x_shift << "\t" << h*v-z_shift; // write the position

            // loop through all OP components...
            for (int n = 0; n < Nop; n++) {
               int id = ID(u, v, n); // get the id
               data << "\t" << real(solution(id)) << "  " << imag(solution(id));
            }

            int id = ID(u, v, 0); // get the id
            // add more here if needed 
            data << "\t" << FEdens(id); 
            data << "\t" << FEbulk(id);
            data << "\t" << FEdens(id) - FEbulk(id); // FE gradient part
            data << "\t" << FEdens(id) - FEref(id); // FE relative to some reference 

            data << endl; // finish the line
		FE_xintegrated += h* (FEdens(id)-FEref(id));
         } // u-loop
	   if(Nu>1) data << "\n"; // for gnuplot to make 3D surfaces we need a break 
         fexint << h*v-z_shift << "\t" << FE_xintegrated << endl;
       } // v-loop
	}
	else std::cout << "Unable to open '" << file_name << "' to write vector to file." << std::endl;
}
