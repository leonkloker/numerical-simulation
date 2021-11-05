#pragma once

#include "array2d.h"
#include <array>

/** A field variable is the discretization of a scalar function f(x) with x in the computational domain. 
 * More specifically, a scalar value is stored at discrete nodes/points. 
 * The nodes are arranged in an equidistant mesh with specified mesh width. 
 */
class FieldVariable : public Array2D
{
  public:
    //! constructor
    FieldVariable(std::array < int, 2 > size, std::array< double, 2 > origin, std::array< double, 2 > meshWidth);

    // Get the value at the Cartesian coordinate (x,y). The value is linearly interpolated between stored points. 
    double interpolateAt (double x, double y) const;

  private:
    const std::array< double, 2 > origin_; // Cartesian coordinates of the point with (i,j) = (0,0), this is different from (0,0) for the u,v and p field variables. 
    const std::array< double, 2 > meshWidth_; // The length of a single cell in x and y direction.
};
