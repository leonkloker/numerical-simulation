#include "cg_parallel.h"

CGParallel::CGParallel(std::shared_ptr<Discretization>discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations) : 
PressureSolver(discretization, partition, epsilon, maximumNumberOfIterations){}

void CGParallel::solve(){
  double dx2 = pow(discretization_->dx(), 2);
  double dy2 = pow(discretization_->dy(), 2);
  int N = partition_.nCellsGlobal()[0] * partition_.nCellsGlobal()[1];

  calculateResidual(dx2, dy2);
  double alpha = 0;
  double beta = 0;
  updateDirection(beta);

  double initial_residual = getGlobalResidual(N);
  double current_residual = initial_residual;
  double old_residual;
  int iteration = 0;

  while (current_residual / initial_residual > epsilon_ && iteration < maximumNumberOfIterations_){
      exchangeDirection();
      updateQ(dx2, dy2);
      alpha = getAlpha(current_residual);
      updatePressure(alpha);
      updateResidual(alpha);
      old_residual = current_residual;
      current_residual = getGlobalResidual(N);
      beta = getBeta(old_residual, current_residual);
      updateDirection(beta);
      setBoundaryValues();
      iteration++;
  }
}

void CGParallel::setBoundaryValues(){
  if (partition_.boundaryTop()){
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++){
      discretization_->p(i, discretization_->pJEnd()-1) = discretization_->p(i, discretization_->pJEnd()-2);
    }
  }
  if (partition_.boundaryBottom()){
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++){
      discretization_->p(i, discretization_->pJBegin()) = discretization_->p(i, discretization_->pJBegin()+1);
    }
  }
  if (partition_.boundaryRight()){
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++){
      discretization_->p(discretization_->pIEnd()-1, j) = discretization_->p(discretization_->pIEnd()-2,j);
    }
  }
  if (partition_.boundaryLeft()){
    for (int j = discretization_->pJBegin(); j < discretization_->pJEnd(); j++){
      discretization_->p(discretization_->pIBegin(),j) = discretization_->p(discretization_->pIBegin()+1,j);
    }
  }
}

double CGParallel::getAlpha(double residual){
  double local_dq = 0;

  for (int i = discretization_->dIBegin() + 1; i < discretization_->dIEnd() - 1; i++){
    for (int j = discretization_->dJBegin() + 1; j < discretization_->dJEnd() - 1; j++){
      local_dq = local_dq + discretization_->d(i,j) * discretization_->q(i,j);
    }
  }

  double global_dq;
  MPI_Allreduce(&local_dq, &global_dq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return pow(residual, 2) / global_dq;
}

double CGParallel::getBeta(double old_residual, double residual){
  return pow(residual / old_residual, 2);
}

double CGParallel::getGlobalResidual(int N){
  double localRes = 0;

  for (int i = discretization_->rIBegin() + 1; i < discretization_->rIEnd() - 1; i++){
    for (int j = discretization_->rJBegin() + 1; j < discretization_->rJEnd() - 1; j++){
      localRes = localRes + pow(discretization_->r(i,j), 2);
    }
  }

  double globalRes;
  MPI_Allreduce(&localRes, &globalRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return pow(globalRes / N, 0.5);
}

void CGParallel::calculateResidual(double dx2, double dy2){
  for (int i = discretization_->rIBegin() + 1; i < discretization_->rIEnd() - 1; i++){
    for (int j = discretization_->rJBegin() + 1; j < discretization_->rJEnd() - 1; j++){
      discretization_->r(i,j) = -((discretization_->p(i-1,j) - 2 * discretization_->p(i,j) + discretization_->p(i+1,j)) / dx2) - 
      ((discretization_->p(i,j-1) - 2 * discretization_->p(i,j) + discretization_->p(i,j+1)) / dy2) + discretization_->rhs(i,j);
    }
  }
}

void CGParallel::updateDirection(double beta){
  for (int i = discretization_->dIBegin() + 1; i < discretization_->dIEnd() - 1; i++){
    for (int j = discretization_->dJBegin() + 1; j < discretization_->dJEnd() - 1; j++){
      discretization_->d(i,j) = discretization_->r(i,j) + beta * discretization_->d(i,j);
    }
  }
}

void CGParallel::updateResidual(double alpha){
  for (int i = discretization_->rIBegin() + 1; i < discretization_->rIEnd() - 1; i++){
    for (int j = discretization_->rJBegin() + 1; j < discretization_->rJEnd() - 1; j++){
      discretization_->r(i,j) = discretization_->r(i,j) - alpha * discretization_->q(i,j);
    }
  }
}

void CGParallel::updateQ(double dx2, double dy2){
  for (int i = discretization_->qIBegin() + 1; i < discretization_->qIEnd() - 1; i++){
    for (int j = discretization_->qJBegin() + 1; j < discretization_->qJEnd() - 1; j++){
      discretization_->q(i,j) = ((discretization_->d(i-1,j) - 2 * discretization_->d(i,j) + discretization_->d(i+1,j)) / dx2) + 
      ((discretization_->d(i,j-1) - 2 * discretization_->d(i,j) + discretization_->d(i,j+1)) / dy2);
    }
  }
}

void CGParallel::updatePressure(double alpha){
  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      discretization_->p(i,j) = discretization_->p(i,j) + alpha * discretization_->d(i,j);
    }
  }
}

void CGParallel::exchangeDirection(){
  MPI_Request sendRequestRight;

  if (!partition_.boundaryRight()){
    std::vector<double> sendBufferRight(partition_.nCells()[1], 0);
        
    for (int j = discretization_->dJBegin()+1; j < discretization_->dJEnd()-1; j++){
      sendBufferRight[j-1] = discretization_->d(discretization_->dIEnd()-2,j);
    }    
    MPI_Isend(sendBufferRight.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &sendRequestRight);
  }

  MPI_Request sendRequestLeft;

  if (!partition_.boundaryLeft()){
    std::vector<double> sendBufferLeft(partition_.nCells()[1], 0);
        
    for (int j = discretization_->dJBegin()+1; j < discretization_->dJEnd()-1; j++){
      sendBufferLeft[j-1] = discretization_->d(discretization_->dIBegin()+1,j);
    }    
    MPI_Isend(sendBufferLeft.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &sendRequestLeft);
  }

  if (!partition_.boundaryRight()){
    std::vector<double> receiveBufferRight(partition_.nCells()[1], 0);
    
    MPI_Request receiveRequestRight;

    MPI_Irecv(receiveBufferRight.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &receiveRequestRight);

    MPI_Wait(&receiveRequestRight, MPI_STATUS_IGNORE);

    for (int j = discretization_->dJBegin()+1; j < discretization_->dJEnd()-1; j++){
      discretization_->d(discretization_->dIEnd()-1,j) = receiveBufferRight[j-1];
    }    
  }

  if (!partition_.boundaryLeft()){
    std::vector<double> receiveBufferLeft(partition_.nCells()[1], 0);
    
    MPI_Request receiveRequestLeft;

    MPI_Irecv(receiveBufferLeft.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &receiveRequestLeft);

    MPI_Wait(&receiveRequestLeft, MPI_STATUS_IGNORE);

    for (int j = discretization_->dJBegin()+1; j < discretization_->dJEnd()-1; j++){
      discretization_->d(discretization_->dIBegin(),j) = receiveBufferLeft[j-1];
    }    
  }

  if (!partition_.boundaryLeft()){
    MPI_Wait(&sendRequestLeft, MPI_STATUS_IGNORE);
  }
  if (!partition_.boundaryRight()){
    MPI_Wait(&sendRequestRight, MPI_STATUS_IGNORE);
  }

  // vertical update //

  MPI_Request sendRequestTop;

  if (!partition_.boundaryTop()){
    std::vector<double> sendBufferTop(partition_.nCells()[0], 0);
        
    for (int i = discretization_->dIBegin()+1; i < discretization_->dIEnd()-1; i++){
      sendBufferTop[i-1] = discretization_->d(i,discretization_->dJEnd()-2);
    }    
    MPI_Isend(sendBufferTop.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &sendRequestTop);
  }

  MPI_Request sendRequestBottom;

  if (!partition_.boundaryBottom()){
    std::vector<double> sendBufferBottom(partition_.nCells()[0], 0);
        
    for (int i = discretization_->dIBegin()+1; i < discretization_->dIEnd()-1; i++){
      sendBufferBottom[i-1] = discretization_->d(i,discretization_->dJBegin()+1);
    }    
    MPI_Isend(sendBufferBottom.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &sendRequestBottom);
  }

  if (!partition_.boundaryTop()){
    std::vector<double> receiveBufferTop(partition_.nCells()[0], 0);
    
    MPI_Request receiveRequestTop;

    MPI_Irecv(receiveBufferTop.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &receiveRequestTop);

    MPI_Wait(&receiveRequestTop, MPI_STATUS_IGNORE);

    for (int i = discretization_->dIBegin()+1; i < discretization_->dIEnd()-1; i++){
      discretization_->d(i,discretization_->dJEnd()-1) = receiveBufferTop[i-1];
    }    
  }

  if (!partition_.boundaryBottom()){
    std::vector<double> receiveBufferBottom(partition_.nCells()[0], 0);
    
    MPI_Request receiveRequestBottom;

    MPI_Irecv(receiveBufferBottom.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &receiveRequestBottom);

    MPI_Wait(&receiveRequestBottom, MPI_STATUS_IGNORE);

    for (int i = discretization_->dIBegin()+1; i < discretization_->dIEnd()-1; i++){
      discretization_->d(i,discretization_->dJBegin()) = receiveBufferBottom[i-1];
    }    
  }

  if (!partition_.boundaryTop()){
    MPI_Wait(&sendRequestTop, MPI_STATUS_IGNORE);
  }
  if (!partition_.boundaryBottom()){
    MPI_Wait(&sendRequestBottom, MPI_STATUS_IGNORE);
  }

}

void CGParallel::exchangePressures(){
  MPI_Request sendRequestRight;

  if (!partition_.boundaryRight()){
    std::vector<double> sendBufferRight(partition_.nCells()[1], 0);
        
    for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
      sendBufferRight[j-1] = discretization_->p(discretization_->pIEnd()-2,j);
    }    
    MPI_Isend(sendBufferRight.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &sendRequestRight);
  }

  MPI_Request sendRequestLeft;

  if (!partition_.boundaryLeft()){
    std::vector<double> sendBufferLeft(partition_.nCells()[1], 0);
        
    for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
      sendBufferLeft[j-1] = discretization_->p(discretization_->pIBegin()+1,j);
    }    
    MPI_Isend(sendBufferLeft.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &sendRequestLeft);
  }

  if (!partition_.boundaryRight()){
    std::vector<double> receiveBufferRight(partition_.nCells()[1], 0);
    
    MPI_Request receiveRequestRight;

    MPI_Irecv(receiveBufferRight.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &receiveRequestRight);

    MPI_Wait(&receiveRequestRight, MPI_STATUS_IGNORE);

    for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
      discretization_->p(discretization_->pIEnd()-1,j) = receiveBufferRight[j-1];
    }    
  }

  if (!partition_.boundaryLeft()){
    std::vector<double> receiveBufferLeft(partition_.nCells()[1], 0);
    
    MPI_Request receiveRequestLeft;

    MPI_Irecv(receiveBufferLeft.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &receiveRequestLeft);

    MPI_Wait(&receiveRequestLeft, MPI_STATUS_IGNORE);

    for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
      discretization_->p(discretization_->pIBegin(),j) = receiveBufferLeft[j-1];
    }    
  }

  if (!partition_.boundaryLeft()){
    MPI_Wait(&sendRequestLeft, MPI_STATUS_IGNORE);
  }
  if (!partition_.boundaryRight()){
    MPI_Wait(&sendRequestRight, MPI_STATUS_IGNORE);
  }

  // vertical update //

  MPI_Request sendRequestTop;

  if (!partition_.boundaryTop()){
    std::vector<double> sendBufferTop(partition_.nCells()[0], 0);
        
    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
      sendBufferTop[i-1] = discretization_->p(i,discretization_->pJEnd()-2);
    }    
    MPI_Isend(sendBufferTop.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &sendRequestTop);
  }

  MPI_Request sendRequestBottom;

  if (!partition_.boundaryBottom()){
    std::vector<double> sendBufferBottom(partition_.nCells()[0], 0);
        
    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
      sendBufferBottom[i-1] = discretization_->p(i,discretization_->pJBegin()+1);
    }    
    MPI_Isend(sendBufferBottom.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &sendRequestBottom);
  }

  if (!partition_.boundaryTop()){
    std::vector<double> receiveBufferTop(partition_.nCells()[0], 0);
    
    MPI_Request receiveRequestTop;

    MPI_Irecv(receiveBufferTop.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &receiveRequestTop);

    MPI_Wait(&receiveRequestTop, MPI_STATUS_IGNORE);

    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
      discretization_->p(i,discretization_->pJEnd()-1) = receiveBufferTop[i-1];
    }    
  }

  if (!partition_.boundaryBottom()){
    std::vector<double> receiveBufferBottom(partition_.nCells()[0], 0);
    
    MPI_Request receiveRequestBottom;

    MPI_Irecv(receiveBufferBottom.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &receiveRequestBottom);

    MPI_Wait(&receiveRequestBottom, MPI_STATUS_IGNORE);

    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
      discretization_->p(i,discretization_->pJBegin()) = receiveBufferBottom[i-1];
    }    
  }

  if (!partition_.boundaryTop()){
    MPI_Wait(&sendRequestTop, MPI_STATUS_IGNORE);
  }
  if (!partition_.boundaryBottom()){
    MPI_Wait(&sendRequestBottom, MPI_STATUS_IGNORE);
  }

}