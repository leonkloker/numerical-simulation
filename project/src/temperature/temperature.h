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
  Temperature(std::shared_ptr<Discretization>discretization, Partitioning partition, double re, double pr, std::array<double,2> meshWidth, std::array<double,16> boundaryConditionsT);

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

  // Reynolds number
  double re_;
  
  // Prandtls number
  double pr_;

  // mesh width of the discretization scheme in both spatial directions
  const std::array<double,2> meshWidth_;

  // Boundary conditions for the temperature, where the first 8 entries treat the dirichlet bc and the rest the neumann bc
  const std::array<double,16> boundaryConditionsT_;

};
