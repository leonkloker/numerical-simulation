#pragma once

#include "discretization/discretization.h"
#include "partitioning/partitioning.h"
#include <mpi.h>
#include <cmath>
#include <memory>


// Successive over-relaxation solver.
class SORParallel
{
public:

  // Constructor
  SORParallel(std::shared_ptr<Discretization> discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations, double omega);

  // Solve the system of the Poisson equation for pressure.
  void solve();

private:

  void setBoundaryValues();

  double getGlobalResidual(double dx2, double dy2);

  double getLocalResidual(double dx2, double dy2);

  void exchangePressures();

  void updateGroupBlack(double dx2, double dy2, double factor);

  void updateGroupRed(double dx2, double dy2, double factor);

  // Performe one Gauss Seidel iteration on the entire pressure field.
  void iterationStep(double dx2, double dy2, double factor);
  
  // Relaxation factor
  double omega_;

  std::shared_ptr<Discretization> discretization_;

  Partitioning partition_;

  double epsilon_;

  int maximumNumberOfIterations_;

};
