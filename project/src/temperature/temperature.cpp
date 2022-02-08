#include "temperature.h"

Temperature::Temperature (std::shared_ptr<Discretization>discretization, Partitioning partition) :
discretization_(discretization), partition_(partition)){}


void ComputationParallel::computeTemperature()
{
    // calculate T at the new time
    for (int i = discretization_->pIBegin() + 1; i < discretization_->pIEnd() - 1; i++){
        for (int j = discretization_->pJBegin() + 1; j < discretization_->pJEnd() - 1; j++){
            discretization_->t(i,j) = discretization_->t(i,j) + dt_ * ((1. / (settings_.re * settings_.pr)) * (discretization_->computeD2tDx2(i,j) + discretization_->computeD2tDy2(i,j)) 
            - discretization_->computeDutDx(i,j) - discretization_->computeDvtDy(i,j));
        }
    }

}