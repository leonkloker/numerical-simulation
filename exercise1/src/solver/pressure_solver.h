#pragma once

#include "discretization/discretization.h"
#include <memory>
#include <cmath>

/* Interface for the pressure solver. 
It computes the pressure field variable such that the continuity equation is fulfilled.
*/
class PressureSolver{
public:

  // Constructor
  PressureSolver(std::shared_ptr<Discretization>discretization, double epsilon, int maximumNumberOfIterations);

  // Solve the Poisson problem for the pressure, using the rhs and p field variables in the staggeredGrid.
  virtual void solve() = 0;

protected:

  // Bound value of the residual
  double epsilon_;

  // Maximum number of iterations
  int maximumNumberOfIterations_;

  // Object holding the needed field variables for rhs and p.
  std::shared_ptr<Discretization> discretization_;

  // Set the boundary values to account for homogenous Neumann boundary conditions, 
  // this has to be called after every iteration.
  void setBoundaryValues();

  // Calculate the residual of the pressure Poiss on equation.
  double getResidual(double dx2, double dy2);

};

