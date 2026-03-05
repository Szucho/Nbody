# Nbody
N-body simulation in c++\\
Half-finished unoptimized project of mine that I wrote under a few hours to demonstrate the leapfrog integrator to one of my friends.  

## builds
contains the nbody.exe file, for now I was lazy to create the CMake file of the project so to compile the project paste the following line into the terminal:
clang++ -g -o builds/nbody.exe sources/nbody.cpp sources/init_cond.cpp -Iheaders
or with your prefered compiler (i.e g++)

## headers
Contains the following files:
 1. **init_cond.h:** header file of init_cond.cpp in sources 
 2. **vec_state.h:** defines the Vec struct with some operator overrides, any Vec behaves as a mathematical vector. State is the state vector, consists of a position and a velocity vector. It is split intentionally as this form seemed more comfortable to me when I wrote the integrator. Deriv is the time derivative of the state vector

## scripts
Reminder that the solver was only tested on binary systems, these scripts were written with binaries in mind
1. **plot.py:** plots the trajectory of the reduced particle, the trajectories of the 2 particles and the energy conservation of the integrator
2. **animate.py:** uses the same logic as plot py to create an animation, I used copilot to write it and then fixed the remaining errors in the code

## Sources
1. **init_cond.cpp:** generates initial conditions for binaries, a c++ implementation of the python code of Professor Michela Mapelli
2. **nbody.cpp:** solves the collisional (no softening, O(n^2)) Newtonian n-body problem currently with a simple leapfrog algorithm. The positions of the bodies are then saved in a txt file. The solver is not limited to binaries but can't handle the collisional UV divergence for now - no regularization.

Future plan: writing the generalized Mikkola algorithmic regularization.
