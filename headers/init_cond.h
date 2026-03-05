/*
 Header file of the binary system initial condition generator script
 Bertalan Szuchovszky 12.26.2025
*/

//if INIT_COND_H is not defined, define it (header guard)
#ifndef INIT_COND_H
//define INIT_COND_H to prevent multiple inclusions
#define INIT_COND_H

#include "vec_state.h" //needed as gen_state type is State
//function declarations (only names and params - header)
/*
Only the function that generates the initial conditions of binaries.  
*/

//function that generates the initial conditions

State gen_state(const double&m1,
                const double&m2,
                const double&a,
                const double&e,
                const double&phase,
                const double&G);
#endif
