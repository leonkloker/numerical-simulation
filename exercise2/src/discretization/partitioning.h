#pragma once

#include <array>
#include <cmath>
#include <iostream>

class Partitioning
{
public:

  //! construct the object with given number of cells in x and y direction
  Partitioning(std::array<int,2> nGlobalCells, int world_rank, int world_size);

  //! get number of cells in each coordinate direction 
  const std::array<int,2> nCells() const;
  
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

private:

  // number of cells in x and y direction in current subdomain
  std::array<int,2> nCells_;
 
  // rank of the neighbouring subdomains going clockwise starting at the top
  // -1 if there is no adjacent subdomain but a global boundary
  std::array<int,4> neighbours_;

  // true if there is a global boundary instead of a subdomain,
  // also going clockwise starting at the top
  std::array<bool,4> boundaries_ = {0};

  // amount of subdomains in the x and y direction
  std::array<int, 2> partitionSize_;

  // rank of the current process
  int rank_;
};
