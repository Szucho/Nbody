#include <cmath>
#include "vec_state.h" //I work with Vec, State structs, these are defined in this header file
#include "init_cond.h" //this cpp file is the source script of the func defined in this header

/* 
 Script to generate initial conditions for binaries
 Bertalan Szuchovszky 12.26.2025.

 The script was written based on the init_cond.py script made by Michela Mapelli
 I decided to create a header file which will allow me to use gen_state in the nbody script
*/



State gen_state(const double&m1,
                const double&m2,
                const double&a,
                const double&e,
                const double&phase,
                const double&G){

  const double mtot = m1 + m2;
  State binary;
  binary.positions.resize(2); //expects positions of 2 bodies
  binary.velocities.resize(2);

  //positions
  //first body
  binary.positions[0].x = -m2/mtot*a*(1.0-e*e)*std::cos(phase)/(1.0+e*std::cos(phase));
  binary.positions[0].y = -m2/mtot*a*(1.0-e*e)*std::sin(phase)/(1.0+e*std::cos(phase));
  binary.positions[0].z = 0.0;
  //second body
  binary.positions[1].x = m1/mtot*a*(1.0-e*e)*std::cos(phase)/(1.0+e*std::cos(phase));
  binary.positions[1].y = m1/mtot*a*(1.0-e*e)*std::sin(phase)/(1.0+e*std::cos(phase));
  binary.positions[1].z = 0.0;

  //velocities
  //first body
  binary.velocities[0].x = -m2/mtot * (e*std::cos(phase)/(1.0+e*std::cos(phase))-1.0)*std::sin(phase)*(1.0+e*std::cos(phase))*std::sqrt(G*mtot/(a*(1.0-e*e)));
  binary.velocities[0].y = -m2/mtot*(e*std::sin(phase)*std::sin(phase)/(1.0 + e*std::cos(phase)) + std::cos(phase))*(1.0 + e*std::cos(phase))*std::sqrt(G*mtot/(a*(1.0-e*e)));
  binary.velocities[0].z = 0.0;
  //second body
  binary.velocities[1].x = m1/mtot * (e*std::cos(phase)/(1.0+e*std::cos(phase))-1.0)*std::sin(phase)*(1.0+e*std::cos(phase))*std::sqrt(G*mtot/(a*(1.0-e*e)));
  binary.velocities[1].y = m1/mtot*(e*std::sin(phase)*std::sin(phase)/(1.0 + e*std::cos(phase)) + std::cos(phase))*(1.0 + e*std::cos(phase))*std::sqrt(G*mtot/(a*(1.0-e*e)));
  binary.velocities[1].z = 0.0;

  return binary;
}
