#include "central_differences.h"

CentralDifferences::CentralDifferences(std::array<int,2> nCells, std::array<double,2> meshWidth) : 
Discretization(nCells, meshWidth){}

double CentralDifferences::computeDu2Dx (int i, int j) const{
  return (pow(u(i+1,j)+u(i,j), 2) - pow(u(i,j)+u(i-1,j), 2)) / (4 * meshWidth_[0]);
}

double CentralDifferences::computeDuvDx (int i, int j) const{
  return (((v(i,j) + v(i+1,j)) * (u(i,j+1) + u(i,j))) - ((v(i,j) + v(i-1,j)) * (u(i-1,j) + u(i-1,j+1)))) / (4 * meshWidth_[0]);
}

double CentralDifferences::computeDuvDy (int i, int j) const{
  return (((u(i,j+1) + u(i,j)) * (v(i+1,j) + v(i,j))) - ((u(i,j) + u(i,j-1)) * (v(i+1,j-1) + v(i,j-1)))) / (4 * meshWidth_[1]);
}

double CentralDifferences::computeDv2Dy (int i, int j) const{
  return (pow(v(i,j+1) + v(i,j), 2) - pow(v(i,j) + v(i,j-1), 2)) / (4 * meshWidth_[1]);
}
