#pragma once

#include <array>
#include <cmath>
#include <iostream>

class Partitioning
{
public:

  //! Partition the domain consisting of nGlobalCells cells in x and y direction into world_size subdomain
  Partitioning(std::array<int,2> nGlobalCells, int world_rank, int world_size);

  //! get number of cells in the subdomain of this process
  std::array<int, 2> nCells() const;

  //! get number of cells in entire domain
  std::array<int, 2> nGlobalCells() const;

  //! get the rank of the neighbouring subdomains
  std::array<int, 4> neighbours() const;

  int neighbourTop() const;
  
  int neighbourRight() const;

  int neighbourBottom() const;

  int neighbourLeft() const;

  bool boundaryTop() const;

  bool boundaryRight() const;

  bool boundaryBottom() const;

  bool boundaryLeft() const;

  //! get the boundary type (false = inner, true = outer / global boundary) of the 4 boundaries of this subdomain
  std::array<bool, 4> boundaries() const;

  //! get the rank of the process which handles this subdomain
  int rank() const;

private:

  // boolean representing the group of the pressure solver the subdomain is contained in 
  // true for start in bottom left cell, false for opposite pattern
  bool group;

  // number of cells in entire domain
  std::array<int,2> nGlobalCells_;

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
