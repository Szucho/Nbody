#include <iostream>
#include <cmath>
#include <vector> //arrays are actually faster, template<size_t>N needed for static creation at compile time
#include <fstream>
#include "../headers/vec_state.h" //Vec, State, Deriv structs are defined here
#include "../headers/init_cond.h" //gen_state function to generate initial conditions is defined here


/*
 Nbody simulation toy code
 Bertalan Szuchovszky 12.23.2025.
 modified on 12.26.2025:
  ->Vec, State and Deriv moved to vec_state.h header file
  ->init_cond implemented and accessed through own header file

 TO DO: -test for n>2
        -Mikkola regularization (first for n=2)
        -other integrators => energy comparison
        -option to give/change output filename (QoL)

 Below I define the time derivative of the state vector of bodies based on Newtonian gravity,
 set up the nbody simulation and ingegrate via leapfrog integrator scheme.
 The output trajectory is saved into a file
*/



Deriv nbody(const std::vector<Vec>& positions,
            const std::vector<Vec>& velocities,
            const std::vector<double>& masses){
  
  /*
  This function calculates the derivs (v_i, a_i) belonging to the 
  state of the i-th body (r_i,v_i) using Newtonian gravity for all i
  */
  const double G = 1.0; //dimensionless

  const int n = masses.size(); //number of bodies based on masses

  Deriv derivative;
  derivative.velocities = velocities; //dr/dt = v
  derivative.accelerations.resize(n, Vec(0,0,0)); //dv/dt / a

  //calculating accelerations, a_i = (r_j-r_i) * G*m_j /|r_j-r_i|^3
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (i==j){
        continue;
      }
      Vec r_ij = positions[j]+(positions[i]*(-1.0)); //this is r_j-r_i = r_j + (-1*r_i), - operator is not overriden
      double dist = r_ij.norm();
      double dist3 = dist*dist*dist;
      derivative.accelerations[i] += r_ij * masses[j]*(G/dist3);
    }
  }
  return derivative;
}


State leapfrog_step(const State& state,
                    const std::vector<double>& masses,
                    const double& h){
  /*
  This function handles the leapfrog steps after the first kick (Euler half-step).
  We assume that the input state at the first iteration contains x_0 and v_{1/2},
  and this half step shift in the velocities remains unchanged during the integration.
  */

  const int n = masses.size();

  State new_state;
  new_state.positions.resize(n);
  new_state.velocities.resize(n);
  
  Deriv deriv = nbody(state.positions, state.velocities, masses); //the derivatives at the current position

  //velocity at half timestep
  std::vector<Vec> v_half(n);
  for (int i=0; i<n; i++){
    v_half[i] = state.velocities[i] + deriv.accelerations[i] * (h/2.0);
  }

  //update positions
  for (int i=0; i<n; i++){
    new_state.positions[i] = state.positions[i] + v_half[i]*h; //x_{n+1} = x_n + h*v_{n+1/2}

  }

  Deriv new_deriv = nbody(new_state.positions, v_half, masses);

  //update velocities based on the half step
  for (int i=0; i<n; i++){
    new_state.velocities[i] = v_half[i] + new_deriv.accelerations[i]*(h/2.0); //v_{n+1} = v_{n+1/2} + a_{n+1}*h/2
  }

  return new_state;
}


State init_leapfrog(const State& state,
                    const std::vector<double>& masses,
                    const double& h){
  
  const int n = masses.size();

  State init_state = state;

  Deriv init_deriv = nbody(state.positions, state.velocities, masses);

  //first v half step is done via Euler step
  for (int i=0; i<n; i++){
    init_state.velocities[i] = state.velocities[i] + init_deriv.accelerations[i]*(h/2.0); //v_{1/2} = v_0 + h/2*a_0
  }
  return init_state; //this now contains x_0 and v_{1/2}
}


std::vector<State> integrate(State init_state,
                            const std::vector<double>& masses,
                            const double& t_start,
                            const double& t_end,
                            const double& h){

  int n_steps = static_cast<int>((t_end-t_start)/h);
  std::vector<State> trajectory;
  //pre-allocating memory for the vector without creating any elements
  trajectory.reserve(n_steps+1); //size is n+1 as we also save t_start

  //initializing half step
  State state = init_leapfrog(init_state, masses, h);
  trajectory.push_back(init_state); //0th step is the initial condition

  //integration loop
  for (int i=0; i<n_steps; i++){
    state = leapfrog_step(state, masses, h); //computing new state via leapfrog step
    trajectory.push_back(state); //storing the step
  }
  return trajectory;
}


int main(){
  
  //init_cond.h includes gen_state which generates initial conditions for binaries based on the following vals
  double m1 = 1.0;
  double m2 = 1.0;
  double a = 1.0;
  double e = 0.5;
  double phi = 0;
  double G = 1.0; //dimensionless

  //type is State which is defined here and in init_cond.cpp aswell, will be cleaned up later on
  State init_cond = gen_state(m1, m2, a, e, phi, G);


  std::cout << "Initial conditions:\n";
  for (int i=0; i<2; i++){
    std::cout << i+1 << "-th body positions\n"
              << init_cond.positions[i].x << " "
              << init_cond.positions[i].y << " "
              << init_cond.positions[i].z << std::endl;
    std::cout << i+1 << "-th body velocities\n"
              << init_cond.velocities[i].x << " "
              << init_cond.velocities[i].y << " "
              << init_cond.velocities[i].z << std::endl;
  };

  //masses
  std::vector<double> masses = {m1, m2};

  //integration parameters
  double t_start =0.0;
  double t_end = 50.0;
  double h = 0.01;
  
  std::vector<State> trajectory = integrate(init_cond, masses, t_start, t_end, h);

  //save trajectory to file
  std::ofstream outfile("trajectory.dat");
  for (size_t i=0; i<masses.size(); i++){
    //first row contains masses
    outfile << masses[i] << " ";
  }
  outfile << "\n";
  for (const auto& state : trajectory) {
    //write positions of both bodies
    outfile << state.positions[0].x << " " 
            << state.positions[0].y << " " 
            << state.positions[0].z << " "
            << state.positions[1].x << " " 
            << state.positions[1].y << " " 
            << state.positions[1].z << " "
            << state.velocities[0].x << " "
            << state.velocities[0].y << " "
            << state.velocities[0].z << " "
            << state.velocities[1].x << " "
            << state.velocities[1].y << " "
            << state.velocities[1].z << std::endl;
  }

  std::cout << "Data saved to trajectory.dat" << std::endl;

  //output to check if everything is aight w/out plotting
  std::cout << "Number of steps: " << trajectory.size() << std::endl;
  std::cout << "Final position of body 1: "
            << trajectory.back().positions[0].x << ", "
            << trajectory.back().positions[0].y << ", "
            << trajectory.back().positions[0].z << std::endl;

  return 0;
}
