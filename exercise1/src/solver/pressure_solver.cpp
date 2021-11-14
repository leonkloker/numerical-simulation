#include "pressure_solver.h"

PressureSolver::PressureSolver(std::shared_ptr<Discretization>discretization, double epsilon, int maximumNumberOfIterations) : 
discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations){}

void PressureSolver::setBoundaryValues(){

  for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++){
    discretization_->p(i,discretization_->pJBegin()) = discretization_->p(i,discretization_->pJBegin()+1);
    discretization_->p(i,discretization_->pJEnd()-1) = discretization_->p(i,discretization_->pJEnd()-2);
  }

  for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++){
    discretization_->p(discretization_->pIBegin(),j) = discretization_->p(discretization_->pIBegin()+1,j);
    discretization_->p(discretization_->pIEnd()-1,j) = discretization_->p(discretization_->pIEnd()-2,j);
  }
}

double PressureSolver::getResidual(double dx2, double dy2){
  double l2 = 0;

  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      l2 = l2 + pow(((discretization_->p(i-1,j) - 2 * discretization_->p(i,j) + discretization_->p(i+1,j)) / dx2) + ((discretization_->p(i,j-1) - 2 * discretization_->p(i,j) + discretization_->p(i,j+1)) / dy2) - discretization_->rhs(i,j), 2);
    }
  }

  return pow(l2, 0.5);
}