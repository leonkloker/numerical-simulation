#pragma once

#include "discretization/discretization.h"
#include <cmath>

class CentralDifferences : public Discretization
{
public:

  //! use the constructor of the base class
  CentralDifferences(std::array<int,2> nCells, std::array<double,2> meshWidth);

  //! compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j) const;

  //! compute the 1st derivative ∂ v^2 / ∂x
  virtual double computeDv2Dy(int i, int j) const;

  //! compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j) const;

  //! compute the 1st derivative ∂ (uv) / ∂y
  virtual double computeDuvDy(int i, int j) const;

  //! compute the 1st derivative ∂ (ut) / ∂x
  virtual double computeDutDx(int i, int j) const;

  //! compute the 1st derivative ∂ (vt) / ∂y
  virtual double computeDvtDy(int i, int j) const;

};
