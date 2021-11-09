#pragma once

#include "discretization/staggered_grid.h"
#include <cmath>

class Discretization :
  public StaggeredGrid
{
public:

  //! construct the object with given number of cells in x and y direction
  Discretization(std::array<int,2> nCells, std::array<double,2> meshWidth);

  //! compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j) const = 0;

  //! compute the 1st derivative ∂ v^2 / ∂x
  virtual double computeDv2Dy(int i, int j) const = 0;

  //! compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j) const = 0;

  //! compute the 1st derivative ∂ (uv) / ∂y
  virtual double computeDuvDy(int i, int j) const = 0;

  //! compute the 2nd derivative ∂ u / ∂x^2
  double computeD2uDx2(int i, int j) const;

  //! compute the 2nd derivative ∂ v / ∂x^2
  double computeD2vDx2(int i, int j) const;

  //! compute the 2nd derivative ∂ u / ∂y^2
  double computeD2uDy2(int i, int j) const;

  //! compute the 2nd derivative ∂ v / ∂y^2
  double computeD2vDy2(int i, int j) const;

  //! compute the 1st derivative ∂ p / ∂x
  double computeDpDx(int i, int j) const;

  //! compute the 1st derivative ∂ p / ∂y
  double computeDpDy(int i, int j) const;

};
