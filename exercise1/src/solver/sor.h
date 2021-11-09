#pragma once

#include "solver/pressure_solver.h"
#include "discretization/discretization.h"
#include <cmath>


class SOR : public PressureSolver{
public:

  SOR(std::shared_ptr<Discretization>discretization, double epsilon, int maximumNumberOfIterations, double omega);

  void solve();

private:

  void iterationStep(double dx2, double dy2, double factor);
  
  double omega_;

};
