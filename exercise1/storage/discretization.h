#pragma once

#include "discretization/0_staggered_grid.h"

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
  virtual double computeDuDx2(int i, int j) const = 0;

  //! compute the 2nd derivative ∂ v / ∂x^2
  virtual double computeDvDx2(int i, int j) const = 0;

  //! compute the 2nd derivative ∂ u / ∂y^2
  virtual double computeDuDy2(int i, int j) const = 0;

  //! compute the 2nd derivative ∂ v / ∂y^2
  virtual double computeDvDy2(int i, int j) const = 0;

}
