#pragma once

#include "solver/pressure_solver_parallel.h"
#include <mpi.h>
#include <cmath>
#include <memory>


//! parallelized successive over-relaxation solver
class CGParallel : public PressureSolver
{
public:

  //! constructor
  CGParallel(std::shared_ptr<Discretization> discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations);

  //! solve the pressure Poisson equation
  void solve() override;

private:

  //! set pressure boundary values 
  void setBoundaryValues();

  //! get the global residual
  double getGlobalResidual(int N);

  //! get alpha for pressure and residual update
  double getAlpha(double residual);

  //! get beta for direction update
  double getBeta(double old_residual, double residual);

  //! calculate the residual of the initial pressure guess
  void calculateResidual(double dx2, double dy2);

  //! update the residual at every grid point in the conjugate gradient iteration
  void updateResidual(double alpha);

  //! update the optimization direction in the conjugate gradient iteration
  void updateDirection(double beta);

  //! update q in the conjugate gradient iteration
  void updateQ(double dx2, double dy2);

  //! update the pressure in the conjugate gradient iteration
  void updatePressure(double alpha);

  //! exchange the values of d at the boundaries of neighbouring subdomains
  void exchangeDirection();

  //! exchange pressure values at the boundaries of neighbouring subdomains
  void exchangePressures();

};
