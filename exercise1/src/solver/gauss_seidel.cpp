#include "solver/gauss_seidel.h"

GaussSeidel::GaussSeidel(std::shared_ptr<Discretization>discretization, double epsilon, int maximumNumberOfIterations) :
PressureSolver(discretization, epsilon, maximumNumberOfIterations){}

void GaussSeidel::solve(){
  double dx2 = pow(discretization_->dx(), 2);
  double dy2 = pow(discretization_->dy(), 2);
  double factor = dx2 * dy2 / (2 * (dx2 + dy2));

  double initial_error = getResidual(dx2, dy2);
  double current_error = initial_error;
  int iteration = 0;

  while (current_error > epsilon_ * initial_error && iteration <= maximumNumberOfIterations_){
      iterationStep(dx2, dy2, factor);
      setBoundaryValues();
      current_error = getResidual(dx2, dy2);
      iteration++;
  }
}

void GaussSeidel::iterationStep(double dx2, double dy2, double factor){

  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      discretization_->p(i,j) = factor * ((discretization_->p(i-1,j) + discretization_->p(i+1,j)) / dx2 + 
      (discretization_->p(i,j-1) + discretization_->p(i,j+1)) / dy2 - discretization_->rhs(i,j));
    }
  }
}
