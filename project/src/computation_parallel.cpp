#include "computation_parallel.h"
#include <mpi.h>
#include <vector>

ComputationParallel::ComputationParallel(){}

void ComputationParallel::initialize(int argc, char* argv[])
{
    // load the settings from the parameter file
    settings_.loadFromFile(argv[1]);

    // calculate dx and dy
    meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];

    // initialize the communicator
    MPI_Init(NULL,NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // partition the global domain into several subdomains
    partition_ = Partitioning(settings_.nCells, world_rank, world_size);

    // set up the discretization scheme
    if (settings_.useDonorCell){
        discretization_ = std::make_shared<DonorCell>(partition_.nCells(), meshWidth_, settings_.alpha);
    }else{
        discretization_ = std::make_shared<CentralDifferences>(partition_.nCells(), meshWidth_);
    }

    // set up the solver for the pressure Poisson equation
    if (settings_.pressureSolver == "SOR"){
        pressureSolver_ = std::make_unique<SORParallel>(discretization_, partition_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    }else{
        pressureSolver_ = std::make_unique<CGParallel>(discretization_, partition_, settings_.epsilon, settings_.maximumNumberOfIterations);       
    }

    // initialize the outputwriters
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partition_);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, partition_);
}

void ComputationParallel::runSimulation()
{
    // initialize time to zero
    double time = 0;

    // write initial state as first output
    applyBoundaryValues();
    outputWriterParaview_->writeFile(time);

    // run the integrator of the Navier-Stokes equations
    while (time + dt_ < settings_.endTime){

        // apply the Dirichlet boundary values
        applyBoundaryValues();

          // calculate the timestep such that the simulation is stable
        computeTimeStepWidth();

        // apply the corresponding boundary values for F and G
        applyBoundaryValuesFG();

        // compute F and G
        computePreliminaryVelocities();

        // exchange F and G at the subdomain boundaries between neighbouring processes
        exchangePreliminaryVelocities();

        // compute the right-hand side of the pressure Poisson equation
        computeRightHandSide();

        // solve the pressure Poisson equation
        computePressure();

        // compute U and V
        computeVelocities();

        // exchange U and V at the subdomain boundaries between neighbouring processes processes
        exchangeVelocities();

        // increase the time
        time = time + dt_;

        // write U, V and P into a vtk file every second
        if (std::fmod(time - dt_, 1) >= std::fmod(time, 1)){
            outputWriterParaview_->writeFile(time);
        }

        // outputWriterText_->writeFile(time);
    }

    //adjust the value of dt such that endTime is reached exactly
    dt_ = settings_.endTime - time;

    // run one more integration step
    if (dt_ > 0.00001){
        applyBoundaryValues();
        computeTimeStepWidth();
        applyBoundaryValuesFG();
        computePreliminaryVelocities();
        exchangePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        exchangeVelocities();

	    time = settings_.endTime;

        outputWriterParaview_->writeFile(time);
    }
}

void ComputationParallel::computeTimeStepWidth()
{
    // temporal stability limit of the diffusion operator
    double diffusionDt = settings_.re * pow(meshWidth_[0] * meshWidth_[1], 2) / (2 * (pow(meshWidth_[0],2) + pow(meshWidth_[1],2)));  

    // temporal stability limit of the convection operator
    double convectionDt = meshWidth_[0] / (discretization_->u().max());
    convectionDt = std::min(convectionDt, meshWidth_[1] / discretization_->v().max());

    // introduce the safety factor tau such that everything is definitely stable
    dt_ = settings_.tau * std::min(convectionDt, diffusionDt);

    // dt is not allowed to exceed maximumDt
    if (settings_.maximumDt < dt_){
        dt_ = settings_.maximumDt;
    }

    // get the minimal dt of all processes
    MPI_Allreduce(&dt_, &dt_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}

void ComputationParallel::applyBoundaryValues()
{   
    // apply the Dirichlet boundary values of U at the bottom and top if there is a global boundary
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++){
        if (partition_.boundaryBottom()){
            discretization_->u(i,discretization_->uJBegin()) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i,discretization_->uJBegin()+1);
        }
        if (partition_.boundaryTop()){
            discretization_->u(i,discretization_->uJEnd()-1) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i,discretization_->uJEnd()-2);
        }
    }  

    // apply the Dirichlet boundary values of V at the bottom and top if there is a global boundary
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++){
        if (partition_.boundaryBottom()){
            discretization_->v(i,discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        }
        if (partition_.boundaryTop()){
            discretization_->v(i,discretization_->vJEnd()-1) = settings_.dirichletBcTop[1];
        }
    }  

    // apply the Dirichlet boundary values of U at the left and right side if there is a global boundary
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++){
        if (partition_.boundaryLeft()){
            discretization_->u(discretization_->uIBegin(),j) = settings_.dirichletBcLeft[0];
        }
        if (partition_.boundaryRight()){
            discretization_->u(discretization_->uIEnd()-1,j) = settings_.dirichletBcRight[0];
        }
    } 
    
    // apply the Dirichlet boundary values of V at the left and right side if there is a global boundary
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++){
        if (partition_.boundaryLeft()){
            discretization_->v(discretization_->vIBegin(),j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin()+1,j);
        }
        if (partition_.boundaryRight()){
            discretization_->v(discretization_->vIEnd()-1,j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd()-2,j);
        }
    }
}

void ComputationParallel::applyBoundaryValuesFG()
{
    // apply the Dirichlet boundary values of F at the bottom and top if there is a global boundary
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++){
        if (partition_.boundaryBottom()){
            discretization_->f(i,discretization_->uJBegin()) = discretization_->u(i,discretization_->uJBegin());
        }
        if (partition_.boundaryTop()){
            discretization_->f(i,discretization_->uJEnd()-1) = discretization_->u(i,discretization_->uJEnd()-1);
        }
    }  

    // apply the Dirichlet boundary values of G at the bottom and top if there is a global boundary
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++){
        if (partition_.boundaryBottom()){
            discretization_->g(i,discretization_->vJBegin()) = discretization_->v(i,discretization_->vJBegin());
        }
        if (partition_.boundaryTop()){
            discretization_->g(i,discretization_->vJEnd()-1) = discretization_->v(i,discretization_->vJEnd()-1);
        }
    }  

    // apply the Dirichlet boundary values of F at the left and right side if there is a global boundary
    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++){
        if (partition_.boundaryLeft()){
            discretization_->f(discretization_->uIBegin(),j) = discretization_->u(discretization_->uIBegin(),j);
        }
        if (partition_.boundaryRight()){     
            discretization_->f(discretization_->uIEnd()-1,j) = discretization_->u(discretization_->uIEnd()-1,j);
        }
    } 
    
    // apply the Dirichlet boundary values of G at the left and right side if there is a global boundary
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++){
        if (partition_.boundaryLeft()){
            discretization_->g(discretization_->vIBegin(),j) = discretization_->v(discretization_->vIBegin(),j);
        }
        if (partition_.boundaryRight()){
            discretization_->g(discretization_->vIEnd()-1,j) = discretization_->v(discretization_->vIEnd()-1,j);
        }
    }
}

void ComputationParallel::computePreliminaryVelocities()
{
    // calculate F in the inner part of the subdomain
    for (int i = discretization_->uIBegin() + 1; i < discretization_->uIEnd() - partition_.boundaryRight(); i++){
        for (int j = discretization_->uJBegin() + 1; j < discretization_->uJEnd() - 1; j++){
            discretization_->f(i,j) = discretization_->u(i,j) + dt_ * ((discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))/settings_.re - discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j) + settings_.g[0]);
        }
    }

    // calculate G in the inner part of the subdomain
    for (int i = discretization_->vIBegin() + 1; i < discretization_->vIEnd() - 1; i++){
        for (int j = discretization_->vJBegin() + 1; j < discretization_->vJEnd() - partition_.boundaryTop(); j++){
            discretization_->g(i,j) = discretization_->v(i,j) + dt_ * ((discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))/settings_.re - discretization_->computeDv2Dy(i,j) - discretization_->computeDuvDx(i,j) + settings_.g[1]);
        }
    }
}

void ComputationParallel::exchangePreliminaryVelocities()
{
    MPI_Request sendRequestTop;

    // send G values at the top boundary to the top neighbour
    if (!partition_.boundaryTop()){
        std::vector<double> sendBufferTop(partition_.nCells()[0], 0);
        
        for (int i = discretization_->vIBegin()+1; i < discretization_->vIEnd()-1; i++){
            sendBufferTop[i-1] = discretization_->g(i,discretization_->vJEnd()-1);
        }    
        MPI_Isend(sendBufferTop.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &sendRequestTop);
    }

    MPI_Request sendRequestRight;

    // send F values at the right boundary to the right neighbour
    if (!partition_.boundaryRight()){
        std::vector<double> sendBufferRight(partition_.nCells()[1], 0);
        
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
            sendBufferRight[j-1] = discretization_->f(discretization_->uIEnd()-1, j);
        }    
        MPI_Isend(sendBufferRight.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &sendRequestRight);
    }

    // receive G values at the bottom boundary from the bottom neighbour
    if (!partition_.boundaryBottom()){
        std::vector<double> receiveBufferBottom(partition_.nCells()[0], 0);

        MPI_Request receiveRequestBottom;
        MPI_Irecv(receiveBufferBottom.data(), partition_.nCells()[0], MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &receiveRequestBottom);
        
        MPI_Wait(&receiveRequestBottom, MPI_STATUS_IGNORE);

        for (int i = discretization_->vIBegin()+1; i < discretization_->vIEnd()-1; i++){
            discretization_->g(i,discretization_->vJBegin()) = receiveBufferBottom[i-1];
        }
    }

    // receive F values at the left boundary from the left neighbour
    if (!partition_.boundaryLeft()){
        std::vector<double> receiveBufferLeft(partition_.nCells()[1], 0);

        MPI_Request receiveRequestLeft;
        MPI_Irecv(receiveBufferLeft.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &receiveRequestLeft);

        MPI_Wait(&receiveRequestLeft, MPI_STATUS_IGNORE);
        
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
            discretization_->f(discretization_->uIBegin(),j) = receiveBufferLeft[j-1];
        }
    }

    if (!partition_.boundaryTop()){
        MPI_Wait(&sendRequestTop, MPI_STATUS_IGNORE);
    }

    if (!partition_.boundaryRight()){
        MPI_Wait(&sendRequestRight, MPI_STATUS_IGNORE);
    }
}

void ComputationParallel::computeRightHandSide()
{
    // calculate the right-hand side in the inner part of the subdomain
    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
        for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
            discretization_->rhs(i,j) = (1 / dt_) * ((discretization_->f(i,j) - discretization_->f(i-1,j)) / discretization_->dx() + (discretization_->g(i,j) - discretization_->g(i,j-1)) / discretization_->dy());
        }
    }
}

void ComputationParallel::computePressure()
{
    // call solver to solve Poisson equation
    pressureSolver_->solve();
}

void ComputationParallel::computeVelocities()
{
    // calculate U at the new time
    for (int i = discretization_->uIBegin()+1; i < discretization_->uIEnd()-partition_.boundaryRight(); i++){
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
                discretization_->u(i,j) = discretization_->f(i,j) - dt_ * discretization_->computeDpDx(i,j);
        }
    }

    //Calculate V at the new time
    for (int i = discretization_->vIBegin()+1; i < discretization_->vIEnd()-1; i++){
        for (int j = discretization_->vJBegin()+1; j < discretization_->vJEnd()-partition_.boundaryTop(); j++){
                discretization_->v(i,j) = discretization_->g(i,j) - dt_ * discretization_->computeDpDy(i,j);
        }
    }
}

void ComputationParallel::exchangeVelocities()
{
    ///////////////////////////////
    /// horizontal value update ///
    ///////////////////////////////

    // send U values to the right
    MPI_Request sendRequestRightU;

    if (!partition_.boundaryRight()){
        std::vector<double> sendBufferRightU(partition_.nCells()[1], 0);
        
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
            sendBufferRightU[j-1] = discretization_->u(discretization_->uIEnd()-1, j);
        }    
        MPI_Isend(sendBufferRightU.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &sendRequestRightU);
    }

    // send U values to the left
    MPI_Request sendRequestLeftU;

    if (!partition_.boundaryLeft()){
        std::vector<double> sendBufferLeftU(partition_.nCells()[1], 0);
        
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
            sendBufferLeftU[j-1] = discretization_->u(discretization_->uIBegin()+1, j);
        }    
        MPI_Isend(sendBufferLeftU.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &sendRequestLeftU);
    }

    // send V values to the right
    MPI_Request sendRequestRightV;

    if (!partition_.boundaryRight()){
        std::vector<double> sendBufferRightV(partition_.nCells()[1], 0);
        
        for (int j = discretization_->vJBegin()+1; j < discretization_->vJEnd(); j++){
            sendBufferRightV[j-1] = discretization_->v(discretization_->vIEnd()-2,j);
        }    
        MPI_Isend(sendBufferRightV.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 1, MPI_COMM_WORLD, &sendRequestRightV);
    }

    // send V values to the left
    MPI_Request sendRequestLeftV;

    if (!partition_.boundaryLeft()){
        std::vector<double> sendBufferLeftV(partition_.nCells()[1], 0);
        
        for (int j = discretization_->vJBegin()+1; j < discretization_->vJEnd(); j++){
            sendBufferLeftV[j-1] = discretization_->v(discretization_->vIBegin()+1,j);
        }    
        MPI_Isend(sendBufferLeftV.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 1, MPI_COMM_WORLD, &sendRequestLeftV);
    }

    // receive U values from the left
    if (!partition_.boundaryLeft()){
        std::vector<double> receiveBufferLeftU(partition_.nCells()[1], 0);

        MPI_Request receiveRequestLeftU;
        MPI_Irecv(receiveBufferLeftU.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 0, MPI_COMM_WORLD, &receiveRequestLeftU);
        
        MPI_Wait(&receiveRequestLeftU, MPI_STATUS_IGNORE);

        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
            discretization_->u(discretization_->uIBegin(),j) = receiveBufferLeftU[j-1];
        }
    }

    // receive U values from the right
    if (!partition_.boundaryRight()){
        std::vector<double> receiveBufferRightU(partition_.nCells()[1], 0);

        MPI_Request receiveRequestRightU;
        MPI_Irecv(receiveBufferRightU.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 0, MPI_COMM_WORLD, &receiveRequestRightU);
        
        MPI_Wait(&receiveRequestRightU, MPI_STATUS_IGNORE);

        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
            discretization_->u(discretization_->uIEnd(),j) = receiveBufferRightU[j-1];
        }
    }

    // receive V values from the left
    if (!partition_.boundaryLeft()){
        std::vector<double> receiveBufferLeftV(partition_.nCells()[1], 0);

        MPI_Request receiveRequestLeftV;
        MPI_Irecv(receiveBufferLeftV.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourLeft(), 1, MPI_COMM_WORLD, &receiveRequestLeftV);
        
        MPI_Wait(&receiveRequestLeftV, MPI_STATUS_IGNORE);

        for (int j = discretization_->vJBegin()+1; j < discretization_->vJEnd(); j++){
            discretization_->v(discretization_->vIBegin(),j) = receiveBufferLeftV[j-1];
        }
    }

    // receive V values from the right
    if (!partition_.boundaryRight()){
        std::vector<double> receiveBufferRightV(partition_.nCells()[1], 0);

        MPI_Request receiveRequestRightV;
        MPI_Irecv(receiveBufferRightV.data(), partition_.nCells()[1], MPI_DOUBLE, partition_.neighbourRight(), 1, MPI_COMM_WORLD, &receiveRequestRightV);
        
        MPI_Wait(&receiveRequestRightV, MPI_STATUS_IGNORE);

        for (int j = discretization_->vJBegin()+1; j < discretization_->vJEnd(); j++){
            discretization_->v(discretization_->vIEnd()-1,j) = receiveBufferRightV[j-1];
        }
    }

    // wait until horizontal exchange of velocities is finished
    if (!partition_.boundaryRight()){
        MPI_Wait(&sendRequestRightU, MPI_STATUS_IGNORE);
        MPI_Wait(&sendRequestRightV, MPI_STATUS_IGNORE);
    }
    if (!partition_.boundaryLeft()){
        MPI_Wait(&sendRequestLeftU, MPI_STATUS_IGNORE);
        MPI_Wait(&sendRequestLeftV, MPI_STATUS_IGNORE);
    }

    /////////////////////////////
    /// vertical value update ///
    /////////////////////////////
    
    // send U values to the top
    MPI_Request sendRequestTopU;

    if (!partition_.boundaryTop()){
        std::vector<double> sendBufferTopU(partition_.nCells()[0]+2, 0);
        
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd()+1; i++){
            sendBufferTopU[i] = discretization_->u(i, discretization_->uJEnd()-2);
        }    
        MPI_Isend(sendBufferTopU.data(), partition_.nCells()[0]+2, MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &sendRequestTopU);
    }

    // send U values to the bottom
    MPI_Request sendRequestBottomU;

    if (!partition_.boundaryBottom()){
        std::vector<double> sendBufferBottomU(partition_.nCells()[0]+2, 0);
        
        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd()+1; i++){
            sendBufferBottomU[i] = discretization_->u(i, discretization_->uJBegin()+1);
        }    
        MPI_Isend(sendBufferBottomU.data(), partition_.nCells()[0]+2, MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &sendRequestBottomU);
    }

    // send V values to the top
    MPI_Request sendRequestTopV;

    if (!partition_.boundaryTop()){
        std::vector<double> sendBufferTopV(partition_.nCells()[0]+2, 0);
        
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++){
            sendBufferTopV[i] = discretization_->v(i, discretization_->vJEnd()-1);
        }    
        MPI_Isend(sendBufferTopV.data(), partition_.nCells()[0]+2, MPI_DOUBLE, partition_.neighbourTop(), 1, MPI_COMM_WORLD, &sendRequestTopV);
    }

    // send V values to the bottom
    MPI_Request sendRequestBottomV;

    if (!partition_.boundaryBottom()){
        std::vector<double> sendBufferBottomV(partition_.nCells()[0]+2, 0);
        
        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++){
            sendBufferBottomV[i] = discretization_->v(i, discretization_->vJBegin()+1);
        }    
        MPI_Isend(sendBufferBottomV.data(), partition_.nCells()[0]+2, MPI_DOUBLE, partition_.neighbourBottom(), 1, MPI_COMM_WORLD, &sendRequestBottomV);
    }

    // receive U values from the top
    if (!partition_.boundaryTop()){
        std::vector<double> receiveBufferTopU(partition_.nCells()[0]+2, 0);

        MPI_Request receiveRequestTopU;
        MPI_Irecv(receiveBufferTopU.data(), partition_.nCells()[0]+2, MPI_DOUBLE, partition_.neighbourTop(), 0, MPI_COMM_WORLD, &receiveRequestTopU);
        
        MPI_Wait(&receiveRequestTopU, MPI_STATUS_IGNORE);

        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd()+1; i++){
            discretization_->u(i, discretization_->uJEnd()-1) = receiveBufferTopU[i];
        }
    }

    // receive U values from the bottom
    if (!partition_.boundaryBottom()){
        std::vector<double> receiveBufferBottomU(partition_.nCells()[0]+2, 0);

        MPI_Request receiveRequestBottomU;
        MPI_Irecv(receiveBufferBottomU.data(), partition_.nCells()[0]+2, MPI_DOUBLE, partition_.neighbourBottom(), 0, MPI_COMM_WORLD, &receiveRequestBottomU);
        
        MPI_Wait(&receiveRequestBottomU, MPI_STATUS_IGNORE);

        for (int i = discretization_->uIBegin(); i < discretization_->uIEnd()+1; i++){
            discretization_->u(i,discretization_->uJBegin()) = receiveBufferBottomU[i];
        }
    }

    // receive V values from the top
    if (!partition_.boundaryTop()){
        std::vector<double> receiveBufferTopV(partition_.nCells()[0]+2, 0);

        MPI_Request receiveRequestTopV;
        MPI_Irecv(receiveBufferTopV.data(), partition_.nCells()[0]+2, MPI_DOUBLE, partition_.neighbourTop(), 1, MPI_COMM_WORLD, &receiveRequestTopV);
        
        MPI_Wait(&receiveRequestTopV, MPI_STATUS_IGNORE);

        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++){
            discretization_->v(i,discretization_->vJEnd()) = receiveBufferTopV[i];
        }
    }

    // receive V values from the bottom
    if (!partition_.boundaryBottom()){
        std::vector<double> receiveBufferBottomV(partition_.nCells()[0]+2, 0);

        MPI_Request receiveRequestBottomV;
        MPI_Irecv(receiveBufferBottomV.data(), partition_.nCells()[0]+2, MPI_DOUBLE, partition_.neighbourBottom(), 1, MPI_COMM_WORLD, &receiveRequestBottomV);
        
        MPI_Wait(&receiveRequestBottomV, MPI_STATUS_IGNORE);

        for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++){
            discretization_->v(i,discretization_->vJBegin()) = receiveBufferBottomV[i];
        }
    }

    // wait until vertical exchange of velocities is finished   
    if (!partition_.boundaryTop()){
        MPI_Wait(&sendRequestTopU, MPI_STATUS_IGNORE);
        MPI_Wait(&sendRequestTopV, MPI_STATUS_IGNORE);
    }
    if (!partition_.boundaryBottom()){
        MPI_Wait(&sendRequestBottomU, MPI_STATUS_IGNORE);
        MPI_Wait(&sendRequestBottomV, MPI_STATUS_IGNORE);
    }
}
