#pragma once

#include "discretization/discretization.h"
#include <memory>

class PressureSolver{
public:

  PressureSolver(std::shared_ptr<Discretization>discretization, double epsilon, int maximumNumberOfIterations);

  virtual void solve();

protected:

  double epsilon_;

  int maximumNumberOfIterations_;

  std::shared_ptr<Discretization> discretization_;

  void setBoundaryValues();

  double getResidual(double dx2, double dy2);

};

