#include "discretization.h"

Discretization::Discretization(std::array<int,2> nCells, std::array<double,2> meshWidth) :
StaggeredGrid(nCells, meshWidth){}

double Discretization::computeD2uDx2(int i, int j) const{
  double d2udx2 = (u(i+1,j) - 2 * u(i,j) + u(i-1,j))/(pow(meshWidth_[0],2));
  return d2udx2;
}

double Discretization::computeD2vDx2(int i, int j) const{
  double d2vdx2 = (v(i+1,j) - 2 * v(i,j) + v(i-1,j))/(pow(meshWidth_[0],2));
  return d2vdx2;
}

double Discretization::computeD2uDy2(int i, int j) const{
  double d2udy2 = (u(i,j+1) - 2 * u(i,j) + u(i,j-1))/(pow(meshWidth_[1],2));
  return d2udy2;
}

double Discretization::computeD2vDy2(int i, int j) const{
  double d2vdy2 = (v(i,j+1) - 2 * v(i,j) + v(i,j-1))/(pow(meshWidth_[1],2));
  return d2vdy2;
}

double Discretization::computeDpDx(int i, int j) const{
  double dpdx = (p(i+1,j) - p(i,j))/meshWidth_[0];
  return dpdx;
}

double Discretization::computeDpDy(int i, int j) const{
  double dpdy = (p(i,j+1) - p(i,j))/meshWidth_[1];
  return dpdy;
}


