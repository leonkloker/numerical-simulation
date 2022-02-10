#include <iostream>
#include <cstdlib>
#include <chrono>

#include "computation_parallel.h"

int main(int argc, char *argv[])
{
  #ifndef NDEBUG
  // only run this code in debug target
  std::cout << "lots of inefficient but informative output . . .";
  #endif

  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }

  ComputationParallel simulation;
  simulation.initialize(argc, argv);
  simulation.runSimulation();

  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
