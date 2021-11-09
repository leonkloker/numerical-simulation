#pragma once

#include "pressure_solver.h"
#include "../discretization/discretization.h"


class GaussSeidel : public PressureSolver{
public:

  GaussSeidel(std::shared_ptr<Discretization>discretization, double epsilon, int maximumNumberOfIterations);

  void solve();

private:

  void iterationStep(double dx2, double dy2, double factor);

};
