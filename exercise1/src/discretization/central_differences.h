#pragma once

#include "discretization/discretization.h"
#include <cmath>

class CentralDifferences : public Discretization
{
public:

  //! use the constructor of the base class
  CentralDifferences(std::array<int,2> nCells, std::array<double,2> meshWidth);

  //! compute the 1st derivative ∂ u^2 / ∂x
  double computeDu2Dx(int i, int j) const override;

  //! compute the 1st derivative ∂ v^2 / ∂x
  double computeDv2Dy(int i, int j) const override;

  //! compute the 1st derivative ∂ (uv) / ∂x
  double computeDuvDx(int i, int j) const override;

  //! compute the 1st derivative ∂ (uv) / ∂y
  double computeDuvDy(int i, int j) const override;

};
