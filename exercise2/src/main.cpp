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

  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;
  auto t1 = high_resolution_clock::now();

  ComputationParallel simulation;
  simulation.initialize(argc, argv);
  simulation.runSimulation();

  auto t2 = high_resolution_clock::now();
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << ms_double.count();

  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
