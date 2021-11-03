#pragma once

#include "discretization/1_discretization.h"

class CentralDifferences : public Discretization
{
public:

  //! use the constructor of the base class
  CentralDifferences(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha);

  //! compute the 1st derivative ∂ u^2 / ∂x
  virtual double computeDu2Dx(int i, int j) const;

  //! compute the 1st derivative ∂ v^2 / ∂x
  virtual double computeDv2Dy(int i, int j) const;

  //! compute the 1st derivative ∂ (uv) / ∂x
  virtual double computeDuvDx(int i, int j) const;

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

  //! compute the 1st derivative ∂ p / ∂x
  virtual double computeDpDx(int i, int j) const = 0;

  //! compute the 1st derivative ∂ p / ∂y
  virtual double computeDpDy(int i, int j) const = 0;
}
