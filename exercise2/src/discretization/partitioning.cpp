#include "partitioning.h"
 
Partitioning::Partitioning(std::array<int,2> nGlobalCells, int world_rank, int world_size){

	// save rank of the current process
	rank_ = world_rank;

	// find the decomposition of the global domain into world_size 
	// subdomains such that the overall boundary between subdomains is minimal
	int minimal_boundary = nGlobalCells[0] * nGlobalCells[1];

	for (double rx = 1; rx <= std::pow(world_size, 0.5); rx++){
		int h = int(rx);
		double ry = world_size / rx;

		if (std::floor(ry) == ry){

			if (nGlobalCells[0] > nGlobalCells[1]){
				rx = ry;
				ry = h;
			}
			int current_boundary = (rx-1)*nGlobalCells[1] + (ry-1)*nGlobalCells[0];

			if (current_boundary < minimal_boundary){
				partitionSize_[0] = rx;
				partitionSize_[1] = ry;
				minimal_boundary = current_boundary;
			}
		}
		rx = h;
	}

	// determine the size of the subdomain such that the overall grid size stays the same
	int leftx = nGlobalCells[0] % partitionSize_[0];
	int lefty = nGlobalCells[1] % partitionSize_[1];

	if (rank_ % partitionSize_[0] <= leftx-1){
		nCells_[0] = nGlobalCells[0] / partitionSize_[0] + 1;
	}else{
		nCells_[0] = nGlobalCells[0] / partitionSize_[0];
	}

	if (std::floor(rank_ / partitionSize_[0]) <= lefty-1){
		nCells_[1] = nGlobalCells[1] / partitionSize_[1] + 1;
	}else{
		nCells_[1] = nGlobalCells[1] / partitionSize_[1];
	}

	// find the rank of the neighbouring subdomains as well as adjacent global boundaries
	neighbours_[0] = rank_ + partitionSize_[0]; 
	neighbours_[1] = rank_ + 1;
	neighbours_[2] = rank_ - partitionSize_[0];
	neighbours_[3] = rank_ - 1; 
 
 	
	if (rank_ % partitionSize_[0] == 0){
		boundaries_[3] = true;
		neighbours_[3] = -1;
	}
	if((rank_ + 1) % partitionSize_[0] == 0){
		boundaries_[1] = true;
		neighbours_[1] = -1;
	}

	if ((std::floor(rank_ / partitionSize_[0]) + 1) == partitionSize_[1]){
		boundaries_[0] = true;
		neighbours_[0] = -1;
	}
	if(rank_ < partitionSize_[0]){
		boundaries_[2] = true;
		neighbours_[2] = -1;
	}
}
