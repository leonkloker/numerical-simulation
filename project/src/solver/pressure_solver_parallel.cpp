#include "pressure_solver_parallel.h"

PressureSolver::PressureSolver(std::shared_ptr<Discretization>discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations) : 
discretization_(discretization), partition_(partition), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations){}