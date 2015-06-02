Project Summery

I have finished all the regular questions of this project and the simulation went well. If it didn't run for the first time, you may try more times and it will show the result. I have Implemented a preconditioner for the conjugate gradient step and output the time for one round of conjugate calculation.

Questions

1. Describe the difference between Lagrangian and Eulerian viewpoints for simulation.  Why is the approach used in this assignment called Semi-Lagragian?

In Lagrangain viewpoint, properties are bounded with particles and particles move around changin their property values. 
In Eulerian viewpoint, properties are bounded with environment, though the position of environment doesn't change, it's properties change.
We are using Eulerian viewpoint in this project but when advecting properties such as velocities, density and temperature, we imagine a particle moving to current location while remaining the original properties, and use it to update the property of current location. As an imaginary particle is brought up, it is called Semi-Lagragian.

2. Smoke and water can both be simulated with fluids.  Briefly explain how they are similar and how they are different. 

Similarities:
They can both be simulated using Eulerian viewpoint. We can split the differential equation and follow the step of update sources, advecting velocity, add external forces, projection etc. 

Differences:
Smoke is under bouyancy force while water isn't.Water is under surface tension while smoke is not.
We need to advect free surface of water.

3. List one advantage and one disadvantage to simulating our fluid on a grid.  Describe two other techniques for simulating fluids and the advantages and disadvantages of each.

Advantage:
It's easy to work with the spatial derivatives such as pressure gradient and viscosity.parts

Disadvantage:
Small features are hard to simulate, they will either arouse sharp high frequency profile or will simple be ignored if we use semi-lagragian.

(i) Smoothed-particle hydrodynamics
This method divides the fluid into a set of discrete particles. These particles have a spatial distance and kernel function is used to smooth the distance. The contributions of each particle to a property are weighted according to their distance from the particle of interest, and their density.

Advantage:
It guarantees conservation of mass without extra computation since the particles themselves represent mass.
Disadvantage:
This method has more difficult rendering process for fluids such as smokes.

(ii) Lattice Boltzmann methods
Assuming fluid consisting of fictive particles, such particles perform consecutive propagation and collision processes over a discrete lattice mesh. 
Advantage:
This method is particle-based and thus are better in dealing with local dynamics such as collision. 
Disadvantage:
Thermo-hydrodynamic scheme is absent
(reference from http://en.wikipedia.org)

4. How do level set surfaces differ from mesh surfaces? List one advantage and disadvantage of using a mesh to represent a surface.  List one advantage and disadvantage to using a Level Set Surface.

Level set surfaces treat the surfaces as a funtion of all zeros and use usigned value to decide whether a position is inside or outside the fluid. Mesh surfaces just assume the fluid is one part and positions outside the surface is not considered in the simulation.

Advantage of Mesh Surface:
Mesh surfaces simpler to simulate.
Disadvantage of Mesh Surface:
The free surfaces cannot be accurately simulated.

Advantage of Level set:
Level set can make operations such as topology changes.
Disadvantage of leval set:
Small features cannot be sampled.