#ifndef VEC_STATE_H
#define VEC_STATE_H

#include <cmath>
#include <vector>

/* 
 Header containing 3d vector algebra, state and deriv structs
 Bertalan Szuchovszky 12.26.2025.

 Vec is a 3D vector struct:
  ->addition, subtraction of 2 vectors: +,-,+=,-= overrides
  ->multiplication of a vector by a const: *scalar
  ->norm of a vector
  ->scalar product: v*u where v,u are Vec types
  ->cross product: v.cross(u) = vxu

 State is a struct that holds the position and velocity vectors
 Deriv is a struct that holds the velocity and acceleration vectors - derivative of state
*/


struct Vec {
  /*
  3d vector structure with operator overrides for vector algebra.
  Basic operations (all we need): addition of 2 vectors, multiplying a vector by a scalar & norm of a vector.
  */
  double x, y, z;

  Vec() : x(0), y(0), z(0) {} //null vector if not specified
  Vec(double x, double y, double z) : x(x), y(y), z(z) {} //x,y,z if specified

  Vec operator + (const Vec& other) const{
    return Vec(x+other.x, y+other.y, z+other.z); 
  }

  Vec operator - (const Vec& other) const{
    return Vec(x-other.x, y-other.y, z-other.z); 
  }

  Vec operator * (double scalar) const{ //multiplying by a scalar
    return Vec(x*scalar, y*scalar, z*scalar);
  }

  //will be needed for the acceleration calculation
  Vec& operator += (const Vec& other) { //addition logic
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

  Vec& operator -= (const Vec& other) { //addition logic
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
 

  double norm() const{ //norm of a vector -> returns a scalar
    return std::sqrt(x*x + y*y + z*z);
  }

  //dot product
  double operator * (const Vec& other) const{
    return x*other.x + y*other.y + z*other.z;
  }

  Vec cross(const Vec& other) const{
    return Vec(y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.z);
  }
};

struct State { //state "vector", contains positions and velocities
  std::vector<Vec> positions;
  std::vector<Vec> velocities;
};

struct Deriv { //time derivative of state vector, contains velocities and accelerations
  std::vector<Vec> velocities;
  std::vector<Vec> accelerations;
};

#endif
