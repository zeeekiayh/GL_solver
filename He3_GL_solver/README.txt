gl_fdm.cpp 
	- main file
he3_classes.cpp
	- this contains code that forms RHS for different OP structures
he3bulk.cpp
	- generic bulk He3 free energy and the derivative w.r.to OP components
	- strong coupling would be implemented here

HEADER files:
readwrite.hpp
	- contains code to read conditions.txt file
structures.hpp
	- where common structures are defined (boundary conditions, conditions_in etc)


-------------------------------------------------------------------------
Original files (not used, but parts are taken into new files)
includes.hpp
	- main GL solver class, 
	- little functions for matrices
	- ID of point between grid and array,
	- boundary cond structure
	- OP class
	- Build_D_Matrices()
	- PlaceSubMatrix with Kronecker tensor product
complex_classes.hpp
	- derived classes for GL solver (3-component, 5 component etc)  
	- the full D-matrix built here
real_classes.hpp
 	- older file (probably unnecessary)
