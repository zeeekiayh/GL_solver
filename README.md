# GL_solver

Using FEM to solve the GL equations.
I am open for any suggestions and help.

Ideas for presentation:

Brainstorming:
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
    - .
    - - What did you do that worked and what did you try that didn't work?
* Conclusions: - What conclusions (if any) did you reach?
  - ...?
* Sources:
  - Anton Vorontsov
  - list all papers
  - cite pictures used
  - Eigen C++ library
