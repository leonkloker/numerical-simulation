#pragma once

#include "discretization/discretization.h"
#include "partitioning/partitioning.h"
#include "solver/pressure_solver_parallel.h"
#include <mpi.h>
#include <cmath>
#include <memory>


//! parallelized successive over-relaxation solver
class SORParallel : public PressureSolver
{
public:

  //! constructor
  SORParallel(std::shared_ptr<Discretization> discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations, double omega);

  //! solve the pressure Poisson equation
  void solve() override;

private:

  //! get the global residual
  double getGlobalResidual(double dx2, double dy2);

  //! get the residual on the subdomain
  double getLocalResidual(double dx2, double dy2);

  //! execute a Gauss-Seidel iteration on one part of the chess pattern
  void updateGroupBlack(double dx2, double dy2, double factor);

  //! execute a Gauss-Seidel iteration on the other part of the chess pattern
  void updateGroupRed(double dx2, double dy2, double factor);

  //! execute a Gauss-Seidel iteration on the entire pressure field.
  void iterationStep(double dx2, double dy2, double factor);
  
  //! relaxation factor
  double omega_;

};
