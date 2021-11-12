#include "donor_cell.h"

DonorCell::DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha) :
Discretization(nCells, meshWidth), alpha_(alpha){}

double DonorCell::computeDu2Dx (int i, int j) const {
  return (1 / (4 * meshWidth_[0])) * (pow(u(i+1,j) + u(i,j), 2) - pow(u(i,j) + u(i-1,j), 2)) + 
  (alpha_ / (4 * meshWidth_[0])) * (abs(u(i+1,j) + u(i,j)) * (u(i,j) - u(i+1,j)) - abs(u(i,j) + u(i-1,j)) * (u(i-1,j) - u(i,j)));
}
double DonorCell::computeDuvDx (int i, int j) const {
  return (1 / (4 * meshWidth_[0])) * (((v(i,j) + v(i+1,j)) * (u(i,j+1) + u(i,j))) - (v(i,j) + v(i-1,j)) * (u(i-1,j) + u(i-1,j+1))) +
  (alpha_ / (4 * meshWidth_[0])) * (abs(u(i,j+1) + u(i,j)) * (v(i,j) - v(i+1,j)) - abs(u(i-1,j) + u(i-1,j+1)) * (v(i-1,j) - v(i,j)));
}

double DonorCell::computeDuvDy (int i, int j) const {
  return (1 / (4 * meshWidth_[1])) * (((v(i,j) + v(i+1,j)) * (u(i,j+1) + u(i,j))) - (v(i,j-1) + v(i+1,j-1)) * (u(i,j) + u(i,j-1))) +
  (alpha_ / (4 * meshWidth_[1])) * (abs(v(i,j) + v(i+1,j)) * (u(i,j) - u(i,j+1)) - abs(v(i,j-1) + v(i+1,j-1)) * (u(i,j-1) - u(i,j)));
}

double DonorCell::computeDv2Dy (int i, int j) const {
  return (1 / (4 * meshWidth_[1])) * (pow(v(i,j+1) + v(i,j), 2) - pow(v(i,j) + v(i,j-1), 2)) + 
  alpha_ * (1 / (4 * meshWidth_[1])) * (abs(v(i,j+1) + v(i,j)) * (v(i,j) - v(i,j+1)) - abs(v(i,j) + v(i,j-1)) * (v(i,j-1) - v(i,j)));
}
