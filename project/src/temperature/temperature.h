#pragma once

#include "discretization/discretization.h"
#include "partitioning/partitioning.h"
#include <memory>
#include <cmath>
#include <mpi.h>

// Class for the calculation of the temperature.
class Temperature
{
public:

  //! constructor
  Temperature(std::shared_ptr<Discretization>discretization, Partitioning partition, double re, double pr, std::array<double,2> meshWidth, std::array<double,2> dirichletBcBottomT, std::array<double,2> dirichletBcTopT, 
  std::array<double,2> dirichletBcLeftT, std::array<double,2> dirichletBcRightT, std::array<double,2> neumannBcBottomT, std::array<double,2> neumannBcTopT, 
  std::array<double,2> neumannBcLeftT, std::array<double,2> neumannBcRightT);

  //! calculate t at the new time
  void computeTemperature(double dt_);

  //! set temperature boundary values 
  void applyBoundaryValues();

  //! exchange temperature values at subdomain boundaries
  void exchangeTemperatures();

private:

  //! partitioning of the global domain
  Partitioning partition_;

  // Object holding the needed field variables for t.
  std::shared_ptr<Discretization> discretization_;

  // Bound value of the residual
  double re_;
  // Bound value of the residual
  double pr_;

  const std::array<double,2> meshWidth_;

  const std::array<double,2> dirichletBcBottomT_;
  const std::array<double,2> dirichletBcTopT_;
  const std::array<double,2> dirichletBcLeftT_;
  const std::array<double,2> dirichletBcRightT_;
  const std::array<double,2> neumannBcBottomT_;
  const std::array<double,2> neumannBcTopT_;
  const std::array<double,2> neumannBcLeftT_;
  const std::array<double,2> neumannBcRightT_;
};
