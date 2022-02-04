#pragma once

#include "storage/fieldvariable.h"
#include <array>

class StaggeredGrid
{
public:

  //! construct the object with given number of cells in x and y direction
  StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth);

  //! get the mesh width, i.e. the length of a single cell in x and y direction 
  const std::array<double,2> meshWidth() const;

  //! get number of cells in each coordinate direction 
  const std::array<int,2> nCells() const;
  
  //! get a reference to field variable u
  const FieldVariable& u() const;

  //! get a reference to field variable u
  const FieldVariable& v() const;
 
  //! get a reference to field variable u
  const FieldVariable& p() const;

  //! access value of u in element (i,j)
  double u(int i, int j) const;

  //! access value of u in element (x,y)
  double& u(int i, int j);

  //! access value of v in element (i,j)
  double v(int i, int j) const;

  //! access value of v in element (x,y)
  double& v(int i, int j);

  //! access value of p in element (i,j)
  double p(int i, int j) const;

  //! access value of p in element (x,y)
  double& p(int i, int j);

  //! access value of r in element (i,j)
  double r(int i, int j) const;

  //! access value of r in element (x,y)
  double& r(int i, int j);

  //! access value of d in element (i,j)
  double d(int i, int j) const;

  //! access value of d in element (x,y)
  double& d(int i, int j);

  //! access value of q in element (i,j)
  double q(int i, int j) const;

  //! access value of q in element (x,y)
  double& q(int i, int j);

  //!access value of rhs in element (i,j)
  double& rhs(int i, int j);

  //! access value of F in element (i,j)
  double& f(int i, int j);

  //! access value of G in element (i,j)
  double& g(int i, int j);

  //! get the mesh width in x-direction, δx
  double dx() const;

  //! get the mesh width in y-direction, δy
  double dy() const;

  //! first valid index for u in x direction
  int	uIBegin() const;

  //! one after last valid index for u in x direction
  int uIEnd() const;

  //! first valid index for u in y direction
  int uJBegin() const; 

  //! one after last valid index for u in y direction
  int uJEnd() const; 	

  //! first valid index for v in x direction
  int vIBegin() const;

  //! one after last valid index for v in x direction
  int	vIEnd() const;

  //! first valid index for v in y direction
  int	vJBegin() const;

  //! one after last valid index for v in y direction
  int vJEnd() const;

  //! first valid index for p in x direction
  int pIBegin() const;

  //! one after last valid index for p in x direction
  int pIEnd() const;

  //! first valid index for p in y direction
  int pJBegin() const;

  //! one after last valid index for p in y direction
  int pJEnd() const;

  //! first valid index for p in x direction
  int rIBegin() const;

  //! one after last valid index for p in x direction
  int rIEnd() const;

  //! first valid index for p in y direction
  int rJBegin() const;

  //! one after last valid index for p in y direction
  int rJEnd() const;

  //! first valid index for p in x direction
  int dIBegin() const;

  //! one after last valid index for p in x direction
  int dIEnd() const;

  //! first valid index for p in y direction
  int dJBegin() const;

  //! one after last valid index for p in y direction
  int dJEnd() const;

  //! first valid index for p in x direction
  int qIBegin() const;

  //! one after last valid index for p in x direction
  int qIEnd() const;

  //! first valid index for p in y direction
  int qJBegin() const;

  //! one after last valid index for p in y direction
  int qJEnd() const;

protected:

  const std::array<int,2> nCells_;
 
  const std::array<double,2> meshWidth_;
  
  FieldVariable u_;
  
  FieldVariable v_;
  
  FieldVariable p_;

  FieldVariable r_;

  FieldVariable d_;

  FieldVariable q_;
  
  FieldVariable rhs_;
  
  FieldVariable f_;
  
  FieldVariable g_;

};
