#pragma once

#include "partitioning/partitioning.h"
#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"
#include "settings/settings.h"
#include "discretization/central_differences.h"
#include "discretization/donor_cell.h"
#include "solver/sor_parallel.h"
#include "solver/cg_parallel.h"
#include "solver/pressure_solver_parallel.h"
#include "temperature/temperature.h"
#include <memory>
#include <cmath>


/** This class handles the main simulation.
 *  It implements the time stepping scheme,
 *   computes all the terms and calls the pressure solver.
 */
class ComputationParallel
{
public:

  ComputationParallel();

  //! initialize the computation object, parse the settings from file that is given as the only command line argument
  void initialize(int argc, char* argv[]);

  //! run the whole simulation until tend
  void runSimulation();

private:

  //! compute the time step width dt from maximum velocities
  void computeTimeStepWidth();

  //! set boundary values of u and v to correct values
  void applyBoundaryValues();

  //! set boundary values of f and g to correct values 
  void applyBoundaryValuesFG();

  //! compute the preliminary velocities, F and G
  void computePreliminaryVelocities();

  //! exchange the values of f and g at the subdomain boundaries
  void exchangePreliminaryVelocities();

  //! exchange the values of u and v that are needed to calculate f and g in the next timestep
  void exchangeVelocities();

  //! compute the right hand side of the Poisson equation for the pressure
  void computeRightHandSide();

  //! compute the new velocities, u,v, from the 
  //! preliminary velocities, F,G and the pressure, p
  void computeVelocities();

  //! struct which contains all the relevant parameters for the simulation
  Settings settings_;

  //! manages the division of the global domain in subdomains and stores 
  //! all the relevant data needed for parallel computations
  Partitioning partition_;

  //! discretization scheme of the Navier Stokes equations 
  std::shared_ptr<Discretization> discretization_;

  //! object which is responsible for solving the pressure Poisson equation via SOR or CG
  std::unique_ptr<PressureSolver> pressureSolver_;

  //! object which is responsible for the calculation of the temperature
  std::unique_ptr<Temperature>temperature_;

  //! object which writes the values of u, v and p in a vtk file for each timestep
  std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaview_;

  //! object which writes the values of u, v, f, g, rhs and p in a text file for each timestep
  std::unique_ptr<OutputWriterTextParallel> outputWriterText_;

  //! mesh width of the discretization scheme in both spatial directions
  std::array<double,2> meshWidth_;

  //! timestep for the time integration
  double dt_;
};

