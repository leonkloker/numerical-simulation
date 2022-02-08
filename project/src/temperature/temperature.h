#pragma once

#include "discretization/discretization.h"
#include "partitioning/partitioning.h"
#include <memory>
#include <cmath>
#include <mpi.h>

// Interface for the calculation of the temperature.
class Temperature
{
public:

  //! constructor
  Temperature(std::shared_ptr<Discretization>discretization, Partitioning partition);

private:

  //! calculate t at the new time
  void computeTemperature();

  //! set temperature boundary values 
  void applyBoundaryValuesT();

  //! exchange temperature values at subdomain boundaries
  void exchangeTemperatures();

  //! partitioning of the global domain
  Partitioning partition_;

  // Object holding the needed field variables for t.
  std::shared_ptr<Discretization> discretization_;

};
