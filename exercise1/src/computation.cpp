#include "computation.h"

void Computation::initialize(int argc, char* argv[])
{
    //Load the seetings from the parameter file
    settings_.loadFromFile(argv[1]);

    //calculate dx and dy
    meshWidth_[0] = settings_.physicalSize[0]/settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1]/settings_.nCells[1];

    //Set up the discretization scheme
    if (settings_.useDonorCell){
        discretization_ = std::make_shared<DonorCell>(DonorCell(settings_.nCells, meshWidth_, settings_.alpha));
    }else{
        discretization_ = std::make_shared<CentralDifferences>(CentralDifferences(settings_.nCells, meshWidth_));
    }

    //Set up the solver for the pressure poisson equation
    if (settings_.pressureSolver == "SOR"){
        pressureSolver_ = std::make_unique<SOR>(SOR(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega));
    }else{
        pressureSolver_ = std::make_unique<GaussSeidel>(GaussSeidel(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations));
    }

    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(OutputWriterParaview(discretization_));
    outputWriterText_ = std::make_unique<OutputWriterText>(OutputWriterText(discretization_));
}

void Computation::runSimulation()
{
    //initialize time and fieldvariables
    double time = 0;
    applyBoundaryValues();
    computeTimeStepWidth();
    outputWriterParaview_->writeFile(time);
    outputWriterText_->writeFile(time);

    //run the integrator of the Navier-Stokes equations
    while (time + dt_ < settings_.endTime){
        applyBoundaryValuesFG();
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        applyBoundaryValues();

        time = time + dt_;
        computeTimeStepWidth();

        //write the results of the current timestep into a file for visualization with the outputWriter_
        outputWriterParaview_->writeFile(time);
        outputWriterText_->writeFile(time);
    }

    //adjust the value of dt such that endTime is reached exactly
    dt_ = settings_.endTime - time;

    if (dt_ > 0.00001){
        applyBoundaryValuesFG();
        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        applyBoundaryValues();

	    time = settings_.endTime;

        outputWriterParaview_->writeFile(time);
    }
}

void Computation::computeTimeStepWidth()
{
    //Stability limit of the diffusion operator
    double diffusionDt = settings_.re * pow(meshWidth_[0] * meshWidth_[1], 2) / (2 * (pow(meshWidth_[0],2) + pow(meshWidth_[1],2)));

    //Stability limit of the convection operator
    double convectionDt = meshWidth_[0]/(discretization_->u().max());
    convectionDt = std::min(convectionDt, meshWidth_[1]/(discretization_->v().max()));

    //Find dt such that everything is stable
    dt_ = settings_.tau * std::min(convectionDt, diffusionDt);

    //dt is not allowed to exceed maximumDt
    if (settings_.maximumDt < dt_){
        dt_ = settings_.maximumDt;
    }
}

void Computation::applyBoundaryValues()
{
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++){
        discretization_->u(i,discretization_->uJBegin()) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i,discretization_->uJBegin()+1);
        discretization_->u(i,discretization_->uJEnd()-1) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i,discretization_->uJEnd()-2);
    }  

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++){
        discretization_->v(i,discretization_->vJBegin()) = settings_.dirichletBcBottom[1];
        discretization_->v(i,discretization_->vJEnd()-1) = settings_.dirichletBcTop[1];
    }  

    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++){
        discretization_->u(discretization_->uIBegin(),j) = settings_.dirichletBcLeft[0];
        discretization_->u(discretization_->uIEnd()-1,j) = settings_.dirichletBcRight[0];
    } 
    
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++){
        discretization_->v(discretization_->vIBegin(),j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(discretization_->vIBegin()+1,j);
        discretization_->v(discretization_->vIEnd()-1,j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(discretization_->vIEnd()-2,j);
    }
}

void Computation::applyBoundaryValuesFG()
{
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++){
        discretization_->f(i,discretization_->uJBegin()) = discretization_->u(i,discretization_->uJBegin());
        discretization_->f(i,discretization_->uJEnd()-1) = discretization_->u(i,discretization_->uJEnd()-1);
    }  

    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++){
        discretization_->g(i,discretization_->vJBegin()) = discretization_->v(i,discretization_->vJBegin());
        discretization_->g(i,discretization_->vJEnd()-1) = discretization_->v(i,discretization_->vJEnd()-1);
    }  

    for (int j = discretization_->uJBegin(); j < discretization_->uJEnd(); j++){
        discretization_->f(discretization_->uIBegin(),j) = discretization_->u(discretization_->uIBegin(),j);
        discretization_->f(discretization_->uIEnd()-1,j) = discretization_->u(discretization_->uIEnd()-1,j);
    } 
    
    for (int j = discretization_->vJBegin(); j < discretization_->vJEnd(); j++){
        discretization_->g(discretization_->vIBegin(),j) = discretization_->v(discretization_->vIBegin(),j);
        discretization_->g(discretization_->vIEnd()-1,j) = discretization_->v(discretization_->vIEnd()-1,j);
    }
}

void Computation::computePreliminaryVelocities()
{
    for (int i = discretization_->uIBegin()+1; i < discretization_->uIEnd()-1; i++){
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
            discretization_->f(i,j) = discretization_->u(i,j) + dt_ * ((discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j))/settings_.re - discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j) + settings_.g[0]);
        }
    }

    for (int i = discretization_->vIBegin()+1; i < discretization_->vIEnd()-1; i++){
        for (int j = discretization_->vJBegin()+1; j < discretization_->vJEnd()-1; j++){
            discretization_->g(i,j) = discretization_->v(i,j) + dt_ * ((discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))/settings_.re - discretization_->computeDv2Dy(i,j) - discretization_->computeDuvDx(i,j) + settings_.g[1]);
        }
    }
}

void Computation::computeRightHandSide()
{
    for (int i = discretization_->pIBegin()+1; i < discretization_->pIEnd()-1; i++){
        for (int j = discretization_->pJBegin()+1; j < discretization_->pJEnd()-1; j++){
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
    for (int i = discretization_->uIBegin()+1; i < discretization_->uIEnd()-1; i++){
        for (int j = discretization_->uJBegin()+1; j < discretization_->uJEnd()-1; j++){
                discretization_->u(i,j) = discretization_->f(i,j) - dt_ * discretization_->computeDpDx(i,j);
        }
    }

    for (int i = discretization_->vIBegin()+1; i < discretization_->vIEnd()-1; i++){
        for (int j = discretization_->vJBegin()+1; j < discretization_->vJEnd()-1; j++){
                discretization_->v(i,j) = discretization_->g(i,j) - dt_ * discretization_->computeDpDy(i,j);
        }
    }
}
