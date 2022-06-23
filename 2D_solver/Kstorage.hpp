#ifndef Kstorage_hpp_
#define Kstorage_hpp_

// we can store full K-matrices here and supply them to the SC_classes which will pick the 
// required components
// the coordinate system can be any set of orthogonal directions (right-handed triad)
// (u,  w,  v)  can be 
// (x,  y,  z), 
// (r, phi, z), etc

#define K1 1.0
#define K2 1.0
#define K3 1.0
#define K23 2.0
#define K123 3.0 

// K-matrices in basis 
//       (Auu, Aww, Avv, Avu, Auv, Auw, Awu, Awv, Avw)
//       (Axx, Ayy, Azz, Azx, Axz, Axy, Ayx, Ayz, Azy)
//       ( 0 ,  1 ,  2  , 3  , 4  , 5  , 6  , 7  , 8  )

const std::vector<std::vector<double>> Kuu = {
	{K123, 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , K1  , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , K1  , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , K123, 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , K1  , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , K1  , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , K123, 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , K1  , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , K1  }};

const std::vector<std::vector<double>> Kww = { 
	{K1  , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , K123, 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , K1  , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , K1  , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , K1  , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , K123, 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , K1  , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , K1  , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , K123}};

const std::vector<std::vector<double>> Kvv = { 
	{K1  , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , K1  , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , K123, 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , K1  , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , K123, 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , K1  , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , K1  , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , K123, 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , K1  }};

const std::vector<std::vector<double>> Kuw = { 
	{0   , 0   , 0   , 0   , 0   , K2  , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , K3  , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , K2  },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{K3  , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , K2  , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , K3  , 0   , 0   , 0   , 0   , 0   }};

const std::vector<std::vector<double>> Kwu = { 
	{0   , 0   , 0   , 0   , 0   , K3  , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , K2  , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , K3  },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{K2  , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , K3  , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , K2  , 0   , 0   , 0   , 0   , 0   }};

const std::vector<std::vector<double>> Kuv = { 
	{0   , 0   , 0   , 0   , K2  , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , K3  , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , K2  , 0   , 0   , 0   , 0   , 0   , 0   },
	{K3  , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , K2  , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , K3  , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   }};

const std::vector<std::vector<double>> Kvu = { 
	{0   , 0   , 0   , 0   , K3  , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , K2  , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , K3  , 0   , 0   , 0   , 0   , 0   , 0   },
	{K2  , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , K3  , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , K2  , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   }};

const std::vector<std::vector<double>> Kvw = { 
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , K3  , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , K2  },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , K2  , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , K3  , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , K2  , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , K3  , 0   , 0   , 0   , 0   , 0   , 0   }};

const std::vector<std::vector<double>> Kwv = { 
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , K2  , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , K3  },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , K3  , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , K2  , 0   , 0   , 0   , 0   },
	{0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , K3  , 0   , 0   , 0   , 0   , 0   , 0   , 0   },
	{0   , 0   , K2  , 0   , 0   , 0   , 0   , 0   , 0   }};

// and if we are working in (x,z) Cartesian plane with no y-derivatives, 
// and OP components  (Axx, Ayy, Ayz) = (eta_m=0,1,2)
// we  fill gradK[m][n](i,j) matrix [i=x,z][j=x,z] for this SC class 
// using Kuu, Kvv, Kuv, Kvu matrices' entries corresponding to these OP components 
// for example (Ayz is component 2 in this example) so gradK[2][2](x,z) = Kuu[7][7];

inline double Kij(int i, int j, int m, int n) {
	if (i == 0) {
		     if (j == 0) return Kuu[m][n];
		else if (j == 1) return Kuw[m][n];
		else if (j == 2) return Kuv[m][n];
	} else if (i == 1) {
		     if (j == 0) return Kwu[m][n];
		else if (j == 1) return Kww[m][n];
		else if (j == 2) return Kwv[m][n];
	} else if (i == 2) {
		     if (j == 0) return Kvu[m][n];
		else if (j == 1) return Kvw[m][n];
		else if (j == 2) return Kvv[m][n];
	}
	return 0.0;
}

#endif
