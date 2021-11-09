#pragma once

#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "settings/settings.h"
#include "discretization/central_differences.h"
#include "discretization/donor_cell.h"
#include "solver/gauss_seidel.h"
#include "solver/sor.h"
#include <memory>

/** This class handles the main simulation.
 *  It implements the time stepping scheme,
 *   computes all the terms and calls the pressure solver.
 */
class Computation
{
public:

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

  //! compute the right hand side of the Poisson equation for the pressure
  void computeRightHandSide();

  //! solve the Poisson equation for the pressure
  void computePressure();

  //! compute the new velocities, u,v, from the 
  //! preliminary velocities, F,G and the pressure, p
  void computeVelocities();

  Settings settings_;

  std::shared_ptr<Discretization> discretization_;

  std::unique_ptr<PressureSolver> pressureSolver_;

  std::unique_ptr<OutputWriterParaview> outputWriterParaview_;

  std::unique_ptr<OutputWriterText> outputWriterText_;

  std::array<double,2> meshWidth_;

  double dt_;
};

