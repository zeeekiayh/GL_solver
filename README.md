# Ginzburg-Landau equation solver

Using FDM to solve the GL equations in Cartesian and Cylindrical coordinates. I am open to any suggestions and help.

#### How to use:

* Must have Makefile, GCC, python3; download the Eigen C++ library and make it accessible to GCC, i.e. in your ```include``` folder.
* Run the code for a system (in Linux terminal):
  - ```$ make``` (if not made already)
  - ```$ ./gl_fdm <conditions_file_name> [c]``` (optional argument ```c``` to specify that the conditions file is for a cylindrical system.)
  - ```$ make plot<#>[c]``` (see more below in **plotting**)
* Procedure for new systems:
  - Write a new conditions file (follow the format in the other conditions files). (Look at how these files are read in: ```readwrite.hpp```.)
    - If you want to solve a system that is not 1-, 3-, or 5-component (or not the same OP structure), you will have to write a new derived class, following th pattern in ```SC_class.hpp```, and add all appropriate code to ```SC_classes_base.cpp``` and ```SC_classes_derived```.
    - Notes about the conditions file:
      - Be careful to not change the format--it is important!
      - In the subsection ```BOUNDARY CONDITIONS```, ```N``` means a Neumann condition, ```D``` means Dirichlet.
      - The first number is the "slip length"; should be small (```~1e-10```) for D-type, and larger for N-type.
      - The second number is the (estimated) OP value on that boundary. (It will be held fixed for D-type, but is only an estimate for N-type.)
      - There are sufficient comments to describe the other things there.
      - To have the solver use normal relaxation methods rather than the acceleration method, set the parameter ```wait``` equal to ```number loops```.
  - Write a new ```initialGuess``` function, in  ```SC_classes_base.cpp``` or ```SC_classes_derived```, with corresponding additions to ```SC_class.hpp```. You may choose to use a pre-exisiting initial guess function.
  - Add the conditional statements to ```gl_fdm.cpp``` to say what derived class and which initial guess function must be used for your conditions file. (Be careful to clean up memory! It can cause very strange errors if you don't!)
  - Run the code as described above.

* The **plotting** (python) script should just take the file that ```WriteAllToFile()``` wrote to. So, you can either use the quick ones through the Makefile: ```make plot5c```, for example, to plot the solution from a 5-component cylindrical system. When the script runs, it will ask if you are plotting for a cylindrical OP--that is just so that it will have the correct plot titles and axes labels. If you want to run the python script yourself, say ```python3 plot_solution.py <solution_file_name>```. You can also use ```solution_comparison.py``` (may still be incomplete) to see the difference between 2 solutions: ```python3 solution_comparison.py <sol_file_name_1> <sol_file_name_2>```. It will ask if each is a cylindrical system, again for labeling. It will give you an error message if the mesh sizes are not the same; and it will warn you if it thinks that the step sizes are inconsistent between the 2.

#### Abstract

The Ginzburg-Landau (GL) theory is often used to determine the stability of superfluid phases. We are interested in some phases of He-3 that that are manifest in confined geometries. We design and implement the finite element method to solve the coupled differential equations derived from the GL theory. We discuss methods for handling boundary conditions, building the matrix equation, and methods of relaxation.

$^3$He is used to better understand the fundamental physics of superconductors and superfluids. The Ginzburg-Landau (GL) theory is often used to determine the stability of superfluid phases. Using weak coupling parameters for a 3-component order parameter (OP), the GL equations are derived and normalized with respect to the bulk value. The method for handling boundary conditions and the process of building the matrix equation are discussed. The finite difference method is used (implemented in C++ code) to solve the differential equations derived here. A brief description is provided on how to use the code. A short discussion is given on the main ideas of the program and how to modify the code for other systems. The program successfully solves for the OP on reasonably large meshes, in roughly 10 seconds.

#### Description from presentaion:
* What did we want to accomplish? Why did we do this research?
  - We wanted to see if we could explain why the "polka-dot" phase appeared (in the experiment by Levitin et. al, 2019) when the stripped phase was predicted (by Vorontsov and Sauls).
* What did I do? I solved the Ginzburg-Landau Equations for Helium-3.
  - Starting with the simplest case (1D, single component), solved and found free-energy of the state.
  - Did the same with a 1D, 2-component OP.
  - Still working on the 2D, 3-component free-energy.
  - Plan to fix/finish the 5-component Cartesian and Cylindrical solvers.
* Background:
  - Superfluids and Order Parameters:
    - Similar to how superconductors undergo phase transitions to the superconducting state, He-3 transitions to a superfluid state. This state, just as in superconductivity, is a new ground state, i.e. energetically more favorable than the "old" ground state.
    - There are several bulk phases of He-3 (image: http://ltl.tkk.fi/images/archive/ab.jpg). "Bulk" refers to a large body of the fluid--large compared to the scale of the coherence length or penetration depth.
    - To describe the phases or states, we use an order parameter (OP). We represent the OP as a 3x3 matrix; row indexes are for spin, column indexes are for orbitals. In the normal state, the OP is simply 0.
    - Many other phases can be obtained as a result of geometric confinement: this affects different spin and orbital indexes.
  - Ginzburg-Landau Equations:
    - The free-energy functional is given by (...in LaTeX file) from the GL theory
    - Minimizing the free-energy functional to get the GL equations (one can use integral methods or the Euler-Lagrange equation)
    - Free-energy calculation--how to do it and what it means
* Solving the GL equations: - what are competing approaches (if you know)?
  - How to obtain equations given the OP
    - Plug OP into general GL equation
    - Calculate the bulk value
    - Normalize the equations
  - Finite Element Method:
    - Building derivative matrices
    - Boundary conditions--ghost points at edges and corners
    - Building the matrix equation
    - The idea of relaxation to solve the system
  - Programming: C++, Eigen, Python
    - Started out by finding the inverse of the matrix, but then realized that sparse matrices in Eigen didn't have an inverse method. We also figured out that there are other solvers in the Eigen library that don't even need the inverse, and are much faster!
    - The normal relaxation method doesn't always converge very well (very dependent on the relaxation parameter). Anton wrote an acceleration method that is more stable and can converge faster.
    - Because of how the relaxation method works, the initial guess can make a huge difference in the speed of convergence, and whether or not it actually converges. We try to build each guess to match the given BC's and to be close to what we expect the solution to be. ...
    - ...also, we could have a smaller relaxation parameter (or let normal relaxation run a bit longer before starting the acceleration). Depending on the BC's (and other parameters, like size), the relaxation parameter would have to be adjusted to converge properly.
    - I hadn't multiplied the solver matrix by h^2 initially, so it had many very small values, which could have been the source of some problems. Now we have multiplied both sides of the matrix equation by h^2 to make a better matrix.
    - When we solved the gl equations with multiple OP components, we had to be careful with the order in building the matrix and the rhs vector. We first had the components cycle (para, perp, para, perp, ...), but with the new D matrices (next point), that wouldn't work--so we organized them by component: (Axx, ..., Ayy, ..., Azz, ...).
    - The first 2 programs had the boundary condition insertions hard-coded into the buildMatrix() methods. For cleaner code and ease of inserting boundary conditions for larger systems, we decided to just build the D matrices, and insert later (reading the BC's from a file).
    - When constructing the final matrix equation, the solver matrix had to be built from several sparse matrices. However, Eigen doesn't support the comma initializer, so I had tried going from sparse to dense and back to sparse (that cause std::bad_alloc errors); I tried keeping the matrices as vectors of triplets, but the BC's were not getting inserted correctly; I was able to use the un-supported Eigen function 'kroneckerProduct()' which works great!
    - Several times, when my code was not working, it was because I was using equations that had not been normalized!
    - I often started coding for the most general case that I could, but then my code would get too large and complicated, with so many errors hidden inside. I would then simplify, find the errors, and then start to generalize it again.
    - To find errors in the code, I would place print statements all over the place; I would comment parts out to see if that section caused the compilation error; I would carefully adjust parameters to see how it affected the output; I would run through the math and physics again to make sure I understood and was using the right equations.
* Further research:
  - fix/finish the 5-component Cartesian solver (done!)
  - Make the cylindrical solver and calculate free-energy (done!)
  - Figure out how to use the general gl equation form in the code (to calculate the rhs and the free-energy) rather than simplifying everything by hand and normalizing. (done!)
* Sources/References:
  - A. B. Vorontsov
  - Levitin et al. "Evidence for a Spatially Modulated Superfluid Phase of 3He under Confinement." Physical Review Letters, Feb 2019.
  - A. B. Vorontsov and J. A. Sauls. "Crystalline Order in Superfluid 3He Films." Physical Review Letters, Jan 2007.
  - J. J. Wiman and J. A. Sauls. "Superfluid Phases of 3He in a Periodic Confined Geometry." J Low Temp Phys, 175:17–30, Oct 2014.
  - J. J. Wiman and J. A. Sauls. "Superfluid Phases of 3He in Nanoscale Channels." Physical Review B, 92:1–13, Oct 2015.
  - Joshua J. Wiman and J. A. Sauls. "Strong-Coupling and the Stripe Phase of 3He." J Low Temp Phys, 184:1054–1070, Jun 2016.
  - Eigen C++ library (https://eigen.tuxfamily.org/index.php?title=Main_Page)
