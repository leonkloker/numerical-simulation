#include "computation.h"
#include <cmath>

void Computation::initialize(int argc, char* argv[])
{
    settings_.loadFromFile(argv[0]);
    meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];
    discretization_ = std::make_shared<Discretization>(Discretization(settings_.nCells, meshWidth_));

    if (settings_.pressureSolver == "SOR"){
            pressureSolver_ = std::make_unique<SOR>(SOR(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega));
    }else{
            pressureSolver_ = std::make_unique<GaussSeidel>(GaussSeidel(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations));
    }

    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(OutputWriterParaview(discretization_));
    outputWriterText_ = std::make_unique<OutputWriter>(OutputWriter(discretization_));
}

void Computation::runSimulation()
{
    computeTimeStepWidth();
    applyBoundaryValues();

    double time = 0;

    for (int i = 0; i <= floor(settings_.endTime/dt_)+1; i++){
        computePreliminaryVelocities();
        applyBoundaryValuesFG();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        applyBoundaryValues();
	time = time + dt_;

        //write the results of the current timestep into a file for visualization with the outputWriter_
    }
}

void Computation::computeTimeStepWidth()
{

    double diffusionDt = settings_.re * pow(meshWidth_[0] * meshWidth_[1], 2) / (pow(meshWidth_[0],2) + pow(meshWidth_[1],2));

    //Assume for the driven cavity that u and v dont exceed the prescribed velocity at the top boundary -> u_max, v_max <= dirichletBcTop[0]
    double convectionDt = meshWidth_[0]/settings_.dirichletBcTop[0];
    convectionDt = std::min(convectionDt, meshWidth_[1]/settings_.dirichletBcTop[0]);

    dt_ = std::min(convectionDt, diffusionDt);

    if (settings_.maximumDt < dt_){
        dt_ = settings_.maximumDt;
    }
}

void Computation::applyBoundaryValues()
{
    for (int i = 0; i <= settings_.nCells[0]; i++){
        discretization_->u(i,0) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i,1);
        discretization_->u(i,settings_.nCells[1]+1) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i,settings_.nCells[1]);
        
        if (i > 0){
            discretization_->v(i,0) = settings_.dirichletBcBottom[1];
            discretization_->v(i,settings_.nCells[1]) = settings_.dirichletBcTop[1];
        }
    }

    for (int j = 0; j <= settings_.nCells[1]; j++){
        discretization_->v(0,j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(1,j);
        discretization_->v(settings_.nCells[1]+1,j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(settings_.nCells[1],j);
        
        if (j > 0){
            discretization_->u(0,j) = settings_.dirichletBcLeft[0];
            discretization_->u(settings_.nCells[0],j) = settings_.dirichletBcRight[0];
        }
    }
}

void Computation::applyBoundaryValuesFG()
{
    for (int i = 0; i <= settings_.nCells[0]; i++){
        discretization_->f(i,0) = discretization_->u(i,0);
        discretization_->f(i,settings_.nCells[1]+1) = discretization_->u(i,settings_.nCells[1]+1);
        
        if (i > 0){
            discretization_->g(i,0) = discretization_->v(i,0);
            discretization_->g(i,settings_.nCells[1]) = discretization_->v(i,settings_.nCells[1]);
        }
    }

    for (int j = 0; j <= settings_.nCells[1]; j++){
        discretization_->g(0,j) = discretization_->v(0,j);
        discretization_->g(settings_.nCells[1]+1,j) = discretization_->v(settings_.nCells[1]+1,j);
        
        if (j > 0){
            discretization_->f(0,j) = discretization_->u(0,j);
            discretization_->f(settings_.nCells[0],j) = discretization_->u(settings_.nCells[0],j);
        }
    }
}

void Computation::computePreliminaryVelocities()
{
    for (int i = 1; i <= settings_.nCells[0]; i++){
        for (int j = 1; j <= settings_.nCells[1]; j++){
            if (i < settings_.nCells[0]){
                discretization_->f(i,j) = discretization_->u(i,j) + dt_ * ((discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))/settings_.re - discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j));
            }
            if (j < settings_.nCells[1]){
                discretization_->g(i,j) = discretization_->u(i,j) + dt_ * ((discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))/settings_.re - discretization_->computeDuvDx(i,j) - discretization_->computeDv2Dy(i,j));
            }
        }
    }

    applyBoundaryValuesFG();
}

void Computation::computeRightHandSide()
{
    for (int i = 1; i <= settings_.nCells[0]; i++){
        for (int j = 1; j <= settings_.nCells[1]; j++){
            discretization_->rhs(i,j) = (1 / dt_) * ((discretization_->f(i,j) - discretization_->f(i-1,j)) / meshWidth_[0] + (discretization_->g(i,j) - discretization_->g(i,j-1)) / meshWidth_[1]);
        }
    }
}

void Computation::computePressure()
{
    pressureSolver_->solve();
}

void Computation::computeVelocities()
{
    for (int i = 1; i <= settings_.nCells[0]; i++){
        for (int j = 1; j <= settings_.nCells[1]; j++){
            if (i < settings_.nCells[0]){
                discretization_->u(i,j) = discretization_->f(i,j) - dt_ * discretization_->computeDpDx(i,j);
            }
            if (j < settings_.nCells[1]){
                discretization_->v(i,j) = discretization_->g(i,j) - dt_ * discretization_->computeDpDy(i,j);
            }
        }
    }
}
