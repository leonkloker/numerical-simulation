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

  //! constructor
  PressureSolver(std::shared_ptr<Discretization>discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations);

  //! solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid.
  virtual void solve() = 0;

protected:

  //! set pressure boundary values 
  void setBoundaryValues();

  //! exchange pressure values at subdomain boundaries
  void exchangePressures();

  //! upper bound for the residual
  double epsilon_;

  //! maximum number of iterations
  int maximumNumberOfIterations_;

  //! partitioning of the global domain
  Partitioning partition_;

  //! object accessing the needed field variables for the pressure computation
  std::shared_ptr<Discretization> discretization_;
};

