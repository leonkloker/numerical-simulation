#include "pressure_solver_parallel.h"

PressureSolver::PressureSolver(std::shared_ptr<Discretization>discretization, Partitioning partition, double epsilon, int maximumNumberOfIterations) : 
discretization_(discretization), partition_(partition), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations){}

void PressureSolver::setBoundaryValues(){
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

void PressureSolver::exchangePressures(){
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