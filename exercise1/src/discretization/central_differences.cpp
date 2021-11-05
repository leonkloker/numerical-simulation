#include "central_differences.h"
#include <cmath>

double CentralDifferences::computeDu2Dx (int i, int j){
  Du2Dx = ((u(i+1,j)+u(i,j))^2 - (u(i,j)+u(i-1,j)))^2 / 4dx;
  return Du2Dx;
}
double CentralDifferences::computeDuvDx (int i, int j){
  DuvDx = ((u(i,j+1)+u(i,j))(v(i-1,j)+v(i,j))-(u(i,j)+u(i,j-1))(v(i+1,j-1)+v(i,j-1)))/4dx;
  return DuvDx;
}
double CentralDifferences::computeDuvDy (int i, int j){
  DuvDy = ((u(i,j+1)+u(i,j))(v(i-1,j)+v(i,j))-(u(i,j)+u(i,j-1))(v(i+1,j-1)+v(i,j-1)))/4dy;
  return DuvDy;
}
double CentralDifferences::computeDv2Dy (int i, int j){
  Dv2Dy = ((v(i+1,j)+v(i+1,j-1))^2 - (v(i,j)+v(i,j-1)))^2 / 4dy;
  return Dv2Dy;
}
