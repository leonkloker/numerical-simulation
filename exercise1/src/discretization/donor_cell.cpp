#include "donor_cell.h"
#include <cmath>

DonorCell::DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha) :
Discretization(nCells, meshWidth), alpha_(alpha){}

double DonorCell::computeDu2Dx (int i, int j) const {
  double du2dx = (1. / (4. * this->meshWidth_[0])) * (std::pow(u(i+1,j) + u(i,j), 2) - std::pow(u(i,j) + u(i-1,j), 2)) + (alpha_ / (4. * this->meshWidth_[0])) * (std::abs(u(i+1,j) + u(i,j)) * (u(i,j) - u(i+1,j)) - std::abs(u(i,j) + u(i-1,j)) * (u(i-1,j) - u(i,j)));
  return du2dx;
}
double DonorCell::computeDuvDx (int i, int j) const {
  double duvdx = (1. / (4. * this->meshWidth_[0])) * (((v(i,j) + v(i+1,j)) * (u(i,j+1) + u(i,j))) - (v(i,j) + v(i-1,j)) * (u(i-1,j) + u(i-1,j+1))) + (alpha_ / (4. * this->meshWidth_[0])) * (std::abs(u(i,j+1) + u(i,j)) * (v(i,j) - v(i+1,j)) - std::abs(u(i-1,j) + u(i-1,j+1)) * (v(i-1,j) - v(i,j)));
  return duvdx;
}

double DonorCell::computeDuvDy (int i, int j) const {
  double duvdy = (1. / (4. * this->meshWidth_[1])) * (((v(i,j) + v(i+1,j)) * (u(i,j+1) + u(i,j))) - (v(i,j-1) + v(i+1,j-1)) * (u(i,j) + u(i,j-1))) + (alpha_ / (4. * this->meshWidth_[1])) * (std::abs(v(i,j) + v(i+1,j)) * (u(i,j) - u(i,j+1)) - std::abs(v(i,j-1) + v(i+1,j-1)) * (u(i,j-1) - u(i,j)));
  return duvdy;
}

double DonorCell::computeDv2Dy (int i, int j) const {
  double dv2dy = (1. / (4. * this->meshWidth_[1])) * (std::pow(v(i,j+1) + v(i,j), 2) - std::pow(v(i,j) + v(i,j-1), 2)) + (alpha_ / (4. * this->meshWidth_[1])) * (std::abs(v(i,j+1) + v(i,j)) * (v(i,j) - v(i,j+1)) - std::abs(v(i,j) + v(i,j-1)) * (v(i,j-1) - v(i,j)));
  return dv2dy;
}
