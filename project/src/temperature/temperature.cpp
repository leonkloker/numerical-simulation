#include "temperature.h"

Temperature::Temperature (std::shared_ptr<Discretization>discretization, Partitioning partition) :
discretization_(discretization), partition_(partition){}


void Temperature::computeTemperature()
{    
    // calculate T at the new time
    for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
            discretization_->t(i,j) = discretization_->t(i,j) + dt_ * ((1. / (settings_.re * settings_.pr)) * (discretization_->computeD2tDx2(i,j) + discretization_->computeD2tDy2(i,j)) 
            - discretization_->computeDutDx(i,j) - discretization_->computeDvtDy(i,j));
        }
    }
    // exchange temperatures between the partitions
    exchangeTemperatures();
}

void Temperature::applyBoundaryValuesT(){
  if (partition_.boundaryTop()){
    if (settings_.dirichletBcTopT[1] == 1){
      for (int i = discretization_->tIBegin(); i < discretization_->tIEnd(); i++){
        discretization_->t(i,discretization_->tJEnd()-1) = 2 * settings_.dirichletBcTopT[0] - discretization_->t(i,discretization_->tJEnd()-2);
      }
    }
    else if (settings_.neumannBcTopT[1] == 1){
      for (int i = discretization_->tIBegin(); i < discretization_->tIEnd(); i++){
        discretization_->t(i,discretization_->tJEnd()-1) = - meshWidth_[1] * settings_.neumannBcTopT[0] + discretization_->t(i,discretization_->tJEnd()-2);
      }
    }
    
  }
  if (partition_.boundaryBottom()){
    if (settings_.dirichletBcBottomT[1] == 1){
      for (int i = discretization_->tIBegin(); i < discretization_->tIEnd(); i++){
        discretization_->t(i,discretization_->tJBegin()) = 2 * settings_.dirichletBcBottomT[0] - discretization_->t(i,discretization_->tJBegin()+1);
      }
    }
    else if (settings_.neumannBcBottomT[1] == 1){
      for (int i = discretization_->tIBegin(); i < discretization_->tIEnd(); i++){
        discretization_->t(i,discretization_->tJBegin()) = - meshWidth_[1] * settings_.neumannBcBottomT[0] + discretization_->t(i,discretization_->tJBegin()+1);
      }
    }
  }
  if (partition_.boundaryRight()){
    if (settings_.dirichletBcRightT[1] == 1){
      for (int j = discretization_->tJBegin(); j < discretization_->tJEnd(); j++){
        discretization_->t(discretization_->tIEnd()-1, j) = 2 * settings_.dirichletBcRightT[0] - discretization_->t(discretization_->tIEnd()-2,j);
      }
    }
    else if (settings_.neumannBcRightT[1] == 1){
      for (int j = discretization_->tJBegin(); j < discretization_->tJEnd(); j++){
        discretization_->t(discretization_->tIEnd()-1, j) = - meshWidth_[0] * settings_.neumannBcRightT[0] + discretization_->t(discretization_->tIEnd()-2,j);
      }
    }
  }
  if (partition_.boundaryLeft()){
    if (settings_.dirichletBcLeftT[1] == 1){
      for (int j = discretization_->tJBegin(); j < discretization_->tJEnd(); j++){
        discretization_->t(discretization_->tIBegin(),j) = 2 * settings_.dirichletBcLeftT[0] - discretization_->t(discretization_->tIBegin()+1,j);
      }
    }
    else if (settings_.neumannBcLeftT[1] == 1){
      for (int j = discretization_->tJBegin(); j < discretization_->tJEnd(); j++){
        discretization_->t(discretization_->tIBegin(),j) = - meshWidth_[0] * settings_.neumannBcLeftT[0] + discretization_->t(discretization_->tIBegin()+1,j);
      }
    }    
  }
}

void Temperature::exchangeTemperatures(){
  MPI_Request sendRequestRight;

  if (!partition_.boundaryRight()){
    std::vector<double> sendBufferRight(partition_.nCells()[1], 0);
        
    for (int j = discretization_->tJBegin()+1; j < discretization_->tJEnd()-1; j++){
      sendBufferRight[j-1] = discretization_->t(discretization_->tIEnd()-2,j);
    }    
    MPI_Isend(sendBufferRight.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &sendRequestRight);
  }

  MPI_Request sendRequestLeft;

  if (!partition_.boundaryLeft()){
    std::vector<double> sendBufferLeft(partition_.nCells()[1], 0);
        
    for (int j = discretization_->tJBegin()+1; j < discretization_->tJEnd()-1; j++){
      sendBufferLeft[j-1] = discretization_->t(discretization_->tIBegin()+1,j);
    }    
    MPI_Isend(sendBufferLeft.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &sendRequestLeft);
  }

  if (!partition_.boundaryRight()){
    std::vector<double> receiveBufferRight(partition_.nCells()[1], 0);
    
    MPI_Request receiveRequestRight;

    MPI_Irecv(receiveBufferRight.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &receiveRequestRight);

    MPI_Wait(&receiveRequestRight, MPI_STATUS_IGNORE);

    for (int j = discretization_->tJBegin()+1; j < discretization_->tJEnd()-1; j++){
      discretization_->t(discretization_->tIEnd()-1,j) = receiveBufferRight[j-1];
    }    
  }

  if (!partition_.boundaryLeft()){
    std::vector<double> receiveBufferLeft(partition_.nCells()[1], 0);
    
    MPI_Request receiveRequestLeft;

    MPI_Irecv(receiveBufferLeft.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &receiveRequestLeft);

    MPI_Wait(&receiveRequestLeft, MPI_STATUS_IGNORE);

    for (int j = discretization_->tJBegin()+1; j < discretization_->tJEnd()-1; j++){
      discretization_->t(discretization_->tIBegin(),j) = receiveBufferLeft[j-1];
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
        
    for (int i = discretization_->tIBegin()+1; i < discretization_->tIEnd()-1; i++){
      sendBufferTop[i-1] = discretization_->t(i,discretization_->tJEnd()-2);
    }    
    MPI_Isend(sendBufferTop.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &sendRequestTop);
  }

  MPI_Request sendRequestBottom;

  if (!partition_.boundaryBottom()){
    std::vector<double> sendBufferBottom(partition_.nCells()[0], 0);
        
    for (int i = discretization_->tIBegin()+1; i < discretization_->tIEnd()-1; i++){
      sendBufferBottom[i-1] = discretization_->t(i,discretization_->tJBegin()+1);
    }    
    MPI_Isend(sendBufferBottom.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &sendRequestBottom);
  }

  if (!partition_.boundaryTop()){
    std::vector<double> receiveBufferTop(partition_.nCells()[0], 0);
    
    MPI_Request receiveRequestTop;

    MPI_Irecv(receiveBufferTop.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &receiveRequestTop);

    MPI_Wait(&receiveRequestTop, MPI_STATUS_IGNORE);

    for (int i = discretization_->tIBegin()+1; i < discretization_->tIEnd()-1; i++){
      discretization_->t(i,discretization_->tJEnd()-1) = receiveBufferTop[i-1];
    }    
  }

  if (!partition_.boundaryBottom()){
    std::vector<double> receiveBufferBottom(partition_.nCells()[0], 0);
    
    MPI_Request receiveRequestBottom;

    MPI_Irecv(receiveBufferBottom.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &receiveRequestBottom);

    MPI_Wait(&receiveRequestBottom, MPI_STATUS_IGNORE);

    for (int i = discretization_->tIBegin()+1; i < discretization_->tIEnd()-1; i++){
      discretization_->t(i,discretization_->tJBegin()) = receiveBufferBottom[i-1];
    }    
  }

  if (!partition_.boundaryTop()){
    MPI_Wait(&sendRequestTop, MPI_STATUS_IGNORE);
  }
  if (!partition_.boundaryBottom()){
    MPI_Wait(&sendRequestBottom, MPI_STATUS_IGNORE);
  }

}