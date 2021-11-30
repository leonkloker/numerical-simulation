#include "storage/fieldvariable.h"
#include <cmath>

FieldVariable::FieldVariable(std::array< int, 2 > size, std::array< double, 2 >	origin,	std::array< double, 2 >	meshWidth): 
Array2D(size), origin_(origin), meshWidth_(meshWidth)	
{
	/**
	 *size      The number of entries in x and y direction.
 	 *origin    Cartesian coordinates of the point with (i,j) = (0,0), this is different from (0,0) for the u,v and p field variables.
	 *meshWidth The length of a single cell in x and y direction. 
	 */
}

void FieldVariable::setToZero(){
	std::fill(data_.begin(), data_.end(), 0);
}

double FieldVariable::max() const {
	
	double max_value = 0;

	for(const double &value : data_){
		if (std::abs(value) > max_value){
			max_value = std::abs(value);
		}
	}
	return max_value;
}

double FieldVariable::interpolateAt(double x, double y)	const
{	
	double x_grid = (x - origin_[0]) / meshWidth_[0];
	double y_grid = (y - origin_[1]) / meshWidth_[1];
	int i_left = floor(x_grid);
	int i_right = i_left + 1;
	int j_bottom = floor(y_grid);
	int j_top = j_bottom + 1;

	if (i_left == (size_[0] - 1)){
		i_right = i_left;
		i_left = i_right - 1;
	}
	if (j_bottom == (size_[1] - 1)){
		j_top = j_bottom;
		j_bottom = j_top - 1;
	}
	
	double value1 = (i_right - x_grid) * (*this)(i_left,j_bottom) + (x_grid - i_left) * (*this)(i_right,j_bottom);
	double value2 = (i_right - x_grid) * (*this)(i_left,j_top) + (x_grid - i_left) * (*this)(i_right,j_top);
	double res = (j_top - y_grid) * value1 + (y_grid - j_bottom)*value2;

	return res;
}

std::vector<double> FieldVariable::data(){
  return data_;
}
