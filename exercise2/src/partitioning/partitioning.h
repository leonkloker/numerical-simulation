#pragma once

#include <array>
#include <cmath>
#include <iostream>

class Partitioning
{
public:

  //! default constructor
  Partitioning();

  //! partition the domain consisting of nGlobalCells cells in x and y direction into world_size subdomains
  Partitioning(std::array<int,2> nGlobalCells, int world_rank, int world_size);

  //! get number of cells in the subdomain of this process
  std::array<int, 2> nCells() const;

  //! get number of cells in entire domain
  std::array<int, 2> nCellsGlobal() const;

  //! get the rank of the neighbouring subdomains
  std::array<int, 4> neighbours() const;

  //! get the rank of the top neighbour
  int neighbourTop() const;
  
  //! get the rank of the right neighbour
  int neighbourRight() const;

  //! get the rank of the bottom neighbour
  int neighbourBottom() const;

  //! get the rank of the left neighbour
  int neighbourLeft() const;

  //! get the boundary type at the top
  bool boundaryTop() const;

  //! get the boundary type on the right
  bool boundaryRight() const;

  //! get the boundary type at the bottom
  bool boundaryBottom() const;

  //! get the boundary type on the left
  bool boundaryLeft() const;

  //! get the starting point for the red Gauss-Seidel update
  bool group() const;

  //! get the boundary type (false = inner, true = outer / global boundary) of the 4 boundaries of this subdomain
  std::array<bool, 4> boundaries() const;

  //! get the amount of cells on the left or below this subdomain
  std::array<int, 2> nodeOffset(); 

  //! get the rank of the process which handles this subdomain
  int rank() const;

private:

  //! boolean representing the group of the pressure solver the subdomain is contained in 
  //! true for start in bottom left cell, false for opposite pattern
  bool group_;

  //! number of cells in entire domain
  std::array<int,2> nCellsGlobal_;

  //! number of cells in x and y direction in current subdomain
  std::array<int,2> nCells_;
 
  //! rank of the neighbouring subdomains going clockwise starting at the top
  //! -1 if there is no adjacent subdomain but a global boundary
  std::array<int,4> neighbours_;

  //! true if there is a global boundary instead of a subdomain,
  //! also going clockwise starting at the top
  std::array<bool,4> boundaries_ = {0};

  //! amount of subdomains in the x and y direction
  std::array<int, 2> partitionSize_;

  //! rank of the current process
  int rank_;

  //! amount of cells on the left or below this subdomain
  std::array<int,2> nodeOffset_;
};
