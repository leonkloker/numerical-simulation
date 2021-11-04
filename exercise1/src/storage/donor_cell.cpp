#include "donor_cell.h"
#include <cmath>

double DonorCell::computeDu2Dx (int i, int j){
  Du2Dx = (1/4dx)(((u(i+1,j)+u(i,j))^2)-(u(i,j)+u(i-1,j))^2) + (1/4dx)(|u(i+1,j)+u(i,j)|*(u(i+1,j)+u(i,j))-|u(i,j)+u(i-1,j)|*(u(i,j)+u(i-1,j)))
}

double DonorCell::computeDuvDx (int i, int j){
  DuvDx = (1/dx) ((v(i,j)+v(i+1,j))/2 * (u(i,j+1)+u(i,j))/2 - (v(i,j-1)+v(i+1,j-1))/2 * (u(i,j)+u(i,j-1))/2) + (1/dx) (|(v(i,j)+v(i+1,j))/2| * (u(i,j+1)+u(i,j))/2 - |(v(i,j-1)+v(i+1,j-1))/2| * (u(i,j)+u(i,j-1))/2)
}

double DonorCell::computeDuvDy (int i, int y){
  DuvDy = (1/dy) ((v(i,j)+v(i+1,j))/2 * (u(i,j+1)+u(i,j))/2 - (v(i,j-1)+v(i+1,j-1))/2 * (u(i,j)+u(i,j-1))/2) + (1/dy) (|(v(i,j)+v(i+1,j))/2| * (u(i,j+1)+u(i,j))/2 - |(v(i,j-1)+v(i+1,j-1))/2| * (u(i,j)+u(i,j-1))/2)
}

double DonorCell::computeDv2Dy (int i, int j){
  Dv2Dy = (1/4dy)(((v(i,j)+v(i,j+1))^2)-(v(i,j)+v(i,j-1))^2) + (1/4dx)(|v(i,j)+v(i,j+1)|*(v(i,j)+v(i,j+1))-|v(i,j)+v(i,j-1)|*(v(i,j)+u(i,j-1)))
}
