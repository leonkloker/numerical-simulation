#pragma once

#include "discretization/discretization.h"
#include <memory>
#include <cmath>

class PressureSolver{
public:

  PressureSolver(std::shared_ptr<Discretization>discretization, double epsilon, int maximumNumberOfIterations);

  virtual void solve() = 0;

protected:

  double epsilon_;

  int maximumNumberOfIterations_;

  std::shared_ptr<Discretization> discretization_;

  void setBoundaryValues();

  double getResidual(double dx2, double dy2);

};

