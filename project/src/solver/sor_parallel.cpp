#include "sor_parallel.h"

SORParallel::SORParallel(std::shared_ptr<Discretization>discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations, double omega) :
PressureSolver(discretization, partition, epsilon, maximumNumberOfIterations), omega_(omega){}


void SORParallel::solve(){
  double dx2 = pow(discretization_->dx(), 2);
  double dy2 = pow(discretization_->dy(), 2);
  double factor = dx2 * dy2 / (2 * (dx2 + dy2));

  double initial_error = getGlobalResidual(dx2, dy2);

  double current_error = initial_error;
  int iteration = 0;
  int N = partition_.nCellsGlobal()[0] * partition_.nCellsGlobal()[1];

  while (current_error/pow(N, 0.5) > epsilon_ && iteration < maximumNumberOfIterations_){
      iterationStep(dx2, dy2, factor);
      current_error = getGlobalResidual(dx2, dy2);
      iteration++;
      setBoundaryValues();
  }
}

double SORParallel::getGlobalResidual(double dx2, double dy2){
  double localRes = getLocalResidual(dx2, dy2);
  double globalRes;
  MPI_Allreduce(&localRes, &globalRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return pow(globalRes, 0.5);
}

double SORParallel::getLocalResidual(double dx2, double dy2){
  double l2 = 0;

  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      l2 = l2 + pow(((discretization_->p(i-1,j) - 2 * discretization_->p(i,j) + discretization_->p(i+1,j)) / dx2) + 
      ((discretization_->p(i,j-1) - 2 * discretization_->p(i,j) + discretization_->p(i,j+1)) / dy2) - discretization_->rhs(i,j), 2);
    }
  }
  return l2;
}

void SORParallel::iterationStep(double dx2, double dy2, double factor){
  updateGroupRed(dx2, dy2, factor);
  exchangePressures();
  updateGroupBlack(dx2, dy2, factor);
  exchangePressures();
}

void SORParallel::updateGroupRed(double dx2, double dy2, double factor){
  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      if ((i + j + 1 - partition_.group()) % 2 == 0){
        discretization_->p(i,j) = (1 - omega_) * discretization_->p(i,j) + omega_ * factor * (((discretization_->p(i-1,j) + discretization_->p(i+1,j)) / dx2) + 
        ((discretization_->p(i,j-1) + discretization_->p(i,j+1)) / dy2) - discretization_->rhs(i,j));
      }
    }
  }
}

void SORParallel::updateGroupBlack(double dx2, double dy2, double factor){
  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      if ((i + j - partition_.group()) % 2 == 0){
        discretization_->p(i,j) = (1 - omega_) * discretization_->p(i,j) + omega_ * factor * (((discretization_->p(i-1,j) + discretization_->p(i+1,j)) / dx2) + 
        ((discretization_->p(i,j-1) + discretization_->p(i,j+1)) / dy2) - discretization_->rhs(i,j));
      }
    }
  }
}
