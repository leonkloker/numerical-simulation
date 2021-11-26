#pragma once

#include "solver/pressure_solver.h"
#include "discretization/discretization.h"
#include <cmath>

// Standard Gauss-Seidel solver.
class GaussSeidel : public PressureSolver{
public:

  // Constructor
  GaussSeidel(std::shared_ptr<Discretization>discretization, double epsilon, int maximumNumberOfIterations);

  // Solve the system of the Poisson equation for pressure.
  void solve() override;

private:

  // Performe one Gauss Seidel iteration on the entire pressure field.
  void iterationStep(double dx2, double dy2, double factor);

};
