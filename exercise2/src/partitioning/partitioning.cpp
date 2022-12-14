#include "partitioning.h"
 
Partitioning::Partitioning(){}

Partitioning::Partitioning(std::array<int,2> nCellsGlobal, int world_rank, int world_size){

	// save rank of the current process
	rank_ = world_rank;
	nCellsGlobal_ = nCellsGlobal;

	// find the decomposition of the global domain into world_size 
	// subdomains such that the overall boundary between subdomains is minimal
	int minimal_boundary = nCellsGlobal[0] * nCellsGlobal[1];

	for (double rx = 1; rx <= std::pow(world_size, 0.5); rx++){
		int h = int(rx);
		double ry = world_size / rx;

		if (std::floor(ry) == ry){

			if (nCellsGlobal[0] > nCellsGlobal[1]){
				rx = ry;
				ry = h;
			}
			int current_boundary = (rx-1)*nCellsGlobal[1] + (ry-1)*nCellsGlobal[0];

			if (current_boundary < minimal_boundary){
				partitionSize_[0] = rx;
				partitionSize_[1] = ry;
				minimal_boundary = current_boundary;
			}
		}
		rx = h;
	}

	// determine the size of the subdomain such that the overall grid size stays the same
	int leftx = nCellsGlobal[0] % partitionSize_[0];
	int lefty = nCellsGlobal[1] % partitionSize_[1];

	if (rank_ % partitionSize_[0] <= leftx-1){
		nCells_[0] = nCellsGlobal[0] / partitionSize_[0] + 1;
	}else{
		nCells_[0] = nCellsGlobal[0] / partitionSize_[0];
	}

	if (std::floor(rank_ / partitionSize_[0]) <= lefty-1){
		nCells_[1] = nCellsGlobal[1] / partitionSize_[1] + 1;
	}else{
		nCells_[1] = nCellsGlobal[1] / partitionSize_[1];
	}

	// Determine the pattern for the red-black pressure solver
	int cells_to_left = 0;
	int cells_to_bottom = 0;

	for (int i = 1; i <= (rank_ % partitionSize_[0]); i++){
		if (i <= leftx){	
			cells_to_left = cells_to_left + (nCellsGlobal[0] / partitionSize_[0]) + 1;
		}else{
			cells_to_left = cells_to_left + (nCellsGlobal[0] / partitionSize_[0]);
		}
	}

	for (int i = 1; i <= std::floor(rank_ / partitionSize_[0]); i++){
		if (i <= lefty){	
			cells_to_bottom = cells_to_bottom + (nCellsGlobal[1] / partitionSize_[1]) + 1;
		}else{
			cells_to_bottom = cells_to_bottom + (nCellsGlobal[1] / partitionSize_[1]);
		}
	}

	nodeOffset_ = {cells_to_left, cells_to_bottom};

	if (cells_to_bottom % 2 == 0){
		if (cells_to_left % 2 == 0){
			group_ = true;
		}else{
			group_ = false;
		}
	}else{
		if (cells_to_left % 2 == 0){
			group_ = false;
		}else{
			group_ = true;
		}
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

std::array<int,2> Partitioning::nCells() const{
	return nCells_;
}

std::array<int, 2> Partitioning::nCellsGlobal() const{
	return nCellsGlobal_;
}

std::array<int,2> Partitioning::nodeOffset(){
	return nodeOffset_;
}

std::array<int,4> Partitioning::neighbours() const{
	return neighbours_;
}

int Partitioning::neighbourTop() const{
	return neighbours_[0];
}
  
int Partitioning::neighbourRight() const{
	return neighbours_[1];
}

int Partitioning::neighbourBottom() const{
	return neighbours_[2];
}

int Partitioning::neighbourLeft() const{
	return neighbours_[3];
}

bool Partitioning::boundaryTop() const{
	return boundaries_[0];
}

bool Partitioning::boundaryRight() const{
	return boundaries_[1];
}

bool Partitioning::boundaryBottom() const{
	return boundaries_[2];
}

bool Partitioning::boundaryLeft() const{
	return boundaries_[3];
}

std::array<bool, 4> Partitioning::boundaries() const{
	return boundaries_;
}

int Partitioning::rank() const{
	return rank_;
}

bool Partitioning::group() const{
	return group_;
}
