#pragma once

#include "discretization/discretization.h"
#include "partitioning/partitioning.h"
#include <memory>
#include <cmath>
#include <mpi.h>

/* Interface for the pressure solver. 
It computes the pressure field variable such that the continuity equation is fulfilled.
*/
class PressureSolver{
public:

  // Constructor
  PressureSolver(std::shared_ptr<Discretization>discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations);

  // Solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid.
  virtual void solve() = 0;

protected:

  //! set pressure boundary values 
  void setBoundaryValues();

  //! exchange pressure values at subdomain boundaries
  void exchangePressures();

  // Bound value of the residual
  double epsilon_;

  // Maximum number of iterations
  int maximumNumberOfIterations_;

  //! partitioning of the global domain
  Partitioning partition_;

  // Object holding the needed field variables for rhs and p.
  std::shared_ptr<Discretization> discretization_;
};

