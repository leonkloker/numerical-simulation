#include "sor_parallel.h"

SORParallel::SORParallel(std::shared_ptr<Discretization>discretization, std::shared_ptr<Partitioning>partition, MPI_COMM_World, double epsilon, int maximumNumberOfIterations, double omega) :
discretization_(discretization), partition_(partition), MPI_COMM_WOLRD_(), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations), omega_(omega){}


void SORParallel::solve(){
  double dx2 = pow(discretization_->dx(), 2);
  double dy2 = pow(discretization_->dy(), 2);
  double factor = dx2 * dy2 / (2 * (dx2 + dy2));

  double initial_error = getGlobalResidual(dx2, dy2);
  double current_error = initial_error;
  int iteration = 0;
  int N = partition_->nGlobalCells()[0] * partition_->nGlobalCells()[1];

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

  MPI_Allreduce(&localRes, &globalRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD_);

  return pow(globalRes, 0.5);
}

double SORParallel::getLocalResidual(double dx2, double dy2){
  double l2 = 0;

  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      l2 = l2 + pow(((discretization_->p(i-1,j) - 2 * discretization_->p(i,j) + discretization_->p(i+1,j)) / dx2) + ((discretization_->p(i,j-1) - 2 * discretization_->p(i,j) + discretization_->p(i,j+1)) / dy2) - discretization_->rhs(i,j), 2);
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

void SORParallel::exchangePressures(){
  MPI_REQUEST sendRequestRight;

  if (!partition_->boundaryRight()){
    std::vector<double> sendBufferRight(partition_->nCells()[1], 0);
        
    for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
      sendBufferRight[j-1] = discretization_->p(discretization_->pIEnd()-2,j);
    }    
    MPI_Isend(sendBufferRight.data(), partition_->nCells()[1], MPI_DOUBLE, partition_->neighbourRight(), 0, MPI_COMM_WORLD, &sendRequestRight);
  }

  MPI_REQUEST sendRequestLeft;

  if (!partition_->boundaryLeft()){
    std::vector<double> sendBufferLeft(partition_->nCells()[1], 0);
        
    for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
      sendBufferLeft[j-1] = discretization_->p(discretization_->pIBegin()+1,j);
    }    
    MPI_Isend(sendBufferLeft.data(), partition_->nCells()[1], MPI_DOUBLE, partition_->neighbourLeft(), 0, MPI_COMM_WORLD, &sendRequestLeft);
  }

  if (!partition_->boundaryRight()){
    std::vector<double> receiveBufferRight(partition_->nCells()[1], 0);
    
    MPI_REQUEST receiveRequestRight;

    MPI_Irecv(receiveBufferRight.data(), partition_->nCells()[1], MPI_DOUBLE, partition_->neighbourRight(), 0, MPI_COMM_WORLD, &receiveRequestRight);

    MPI_WAIT(&receiveRequestRight, MPI_STATUS_IGNORE);

    for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
      discretization_->p(discretization_->pIEnd()-1,j) = receiveBufferRight[j-1];
    }    
  }

  if (!partition_->boundaryLeft()){
    std::vector<double> receiveBufferLeft(partition_->nCells()[1], 0);
    
    MPI_REQUEST receiveRequestLeft;

    MPI_Irecv(receiveBufferLeft.data(), partition_->nCells()[1], MPI_DOUBLE, partition_->neighbourLeft(), 0, MPI_COMM_WORLD, &receiveRequestLeft);

    MPI_WAIT(&receiveRequestLeft, MPI_STATUS_IGNORE);

    for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
      discretization_->p(discretization_->pIBegin(),j) = receiveBufferLeft[j-1];
    }    
  }

  MPI_WAIT(&sendRequestLeft, MPI_STATUS_IGNORE);
  MPI_WAIT(&sendRequestRight, MPI_STATUS_IGNORE);

  // vertical update //

  MPI_REQUEST sendRequestTop;

  if (!partition_->boundaryTop()){
    std::vector<double> sendBufferTop(partition_->nCells()[0], 0);
        
    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
      sendBufferRight[i-1] = discretization_->p(i,discretization_->pJEnd()-2);
    }    
    MPI_Isend(sendBufferTop.data(), partition_->nCells()[0], MPI_DOUBLE, partition_->neighbourTop(), 0, MPI_COMM_WORLD, &sendRequestTop);
  }

  MPI_REQUEST sendRequestBottom;

  if (!partition_->boundaryBottom()){
    std::vector<double> sendBufferBottom(partition_->nCells()[0], 0);
        
    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
      sendBufferBottom[i-1] = discretization_->p(i,discretization_->pJBegin()+1);
    }    
    MPI_Isend(sendBufferBottom.data(), partition_->nCells()[0], MPI_DOUBLE, partition_->neighbourBottom(), 0, MPI_COMM_WORLD, &sendRequestBottom);
  }

  if (!partition_->boundaryTop()){
    std::vector<double> receiveBufferTop(partition_->nCells()[0], 0);
    
    MPI_REQUEST receiveRequestTop;

    MPI_Irecv(receiveBufferTop.data(), partition_->nCells()[0], MPI_DOUBLE, partition_->neighbourTop(), 0, MPI_COMM_WORLD, &receiveRequestTop);

    MPI_WAIT(&receiveRequestTop, MPI_STATUS_IGNORE);

    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
      discretization_->p(i,discretization_->pJEnd()-1) = receiveBufferTop[i-1];
    }    
  }

  if (!partition_->boundaryBottom()){
    std::vector<double> receiveBufferBottm(partition_->nCells()[0], 0);
    
    MPI_REQUEST receiveRequestBottom;

    MPI_Irecv(receiveBufferBottom.data(), partition_->nCells()[0], MPI_DOUBLE, partition_->neighbourBottom(), 0, MPI_COMM_WORLD, &receiveRequestBottom);

    MPI_WAIT(&receiveRequestBottom, MPI_STATUS_IGNORE);

    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
      discretization_->p(i,discretization_->pJBegin()) = receiveBufferBottom[i-1];
    }    
  }

  MPI_WAIT(&sendRequestTop, MPI_STATUS_IGNORE);
  MPI_WAIT(&sendRequestBottom, MPI_STATUS_IGNORE);

}

void SORParallel::updateGroupRed(double dx2, double dy2, double factor){
  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      if ((i + j + 1 - partition_->group) % 2 == 0){
        discretization_->p(i,j) = (1 - omega_) * discretization_->p(i,j) + omega_ * factor * (((discretization_->p(i-1,j) + discretization_->p(i+1,j)) / dx2) + 
        ((discretization_->p(i,j-1) + discretization_->p(i,j+1)) / dy2) - discretization_->rhs(i,j));
      }
    }
  }
}

void SORParallel::updateGroupBlack(double dx2, double dy2, double factor){
  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      if ((i + j - partition_->group) % 2 == 0){
        discretization_->p(i,j) = (1 - omega_) * discretization_->p(i,j) + omega_ * factor * (((discretization_->p(i-1,j) + discretization_->p(i+1,j)) / dx2) + 
        ((discretization_->p(i,j-1) + discretization_->p(i,j+1)) / dy2) - discretization_->rhs(i,j));
      }
    }
  }
}


  for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
    for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
      discretization_->p(i,j) = (1 - omega_) * discretization_->p(i,j) + omega_ * factor * (((discretization_->p(i-1,j) + discretization_->p(i+1,j)) / dx2) + 
      ((discretization_->p(i,j-1) + discretization_->p(i,j+1)) / dy2) - discretization_->rhs(i,j));
    }
  }
}
