#include "discretization/staggered_grid.h"
#include "storage/fieldvariable.h"
 
StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth): 
nCells_(nCells),meshWidth_(meshWidth)	
{
	/**
	 *nCells	number of elements in x and y direction 
	 *meshWidth The length of a single cell in x and y direction. 
	 */
}

double StaggeredGrid::dx ()	const {
	// get the mesh width in y-direction, δy 
	return meshWidth_[0];
}

double StaggeredGrid::dy ()	const {
	// get the mesh width in y-direction, δy 
	return meshWidth_[1];
}

double & StaggeredGrid::f(int i, int j){
	// access value of F in element (i,j) 

	// access value of F in element (x,y) 
	return interpolateAt(x,y);
} 		

double & StaggeredGrid::g(int i, int j){
	//Access value of G in element (i,j) 

	// Access value of G in element (x,y) 
	return interpolateAt(x,y);
} 		

const std::array< double, 2 > StaggeredGrid::meshWidth() const {
	// Get the mesh width, i.e. the length of a single cell in x 
	// and y direction 
	return meshWidth_;
}

const std::array< int, 2 > StaggeredGrid::nCells() const {
	// Get number of cells in each coordinate direction 
	return nCells_;
}

const FieldVariable & StaggeredGrid::p() const {
	// Part 1/3
	// Get a reference to field variable u 
}
double StaggeredGrid::p(int	i,int j) const {
	// Part 2/3
	// Access value of p in element (i,j)

	// Get const value of p in element (x,y)

}

double & StaggeredGrid::p(int i, int j) {
	// Part 3/3
	// Access value of p in element (x,y) 
}

int StaggeredGrid::pIBegin() const {
	// First valid index for p in x direction 
}

int StaggeredGrid::pIEnd 	( 		) 	const {
	// One after last valid index for p in x direction 
}

int StaggeredGrid::pJBegin() const {
	// First valid index for p in y direction 
}

int StaggeredGrid::pJEnd() const {
	// One after last valid index for p in y direction 
}

double & StaggeredGrid::rhs(int	i, int j) {
	// Access value of rhs in element (i,j)

	// Access value of p in element (x,y)

}

const FieldVariable & StaggeredGrid::u () const{
	// Part 1/3
	// Get a reference to field variable u 
}

double StaggeredGrid::u (int i,	int	j) const {
	// Part 2/3
	// Access value of u in element (i,j)

	// Get const value of u in element (x,y)

}

double & StaggeredGrid::u (int i, int j) {
	// Access value of u in element (x,y) 
} 		

int StaggeredGrid::uIBegin () const {
	// First valid index for u in x direction 
}

int StaggeredGrid::uIEnd ()	const {
	// One after last valid index for u in x direction 
}

int StaggeredGrid::uJBegin () const {
	// First valid index for u in y direction 
}

int StaggeredGrid::uJEnd () const {
	// One after last valid index for u in y direction 
}

const FieldVariable & StaggeredGrid::v () const {
	// Part 1/3
	// Get a reference to field variable u 
}

double StaggeredGrid::v (int i,	int	j) const {
	// Part 2/3
	// Access value of v in element (i,j)

	// Get const value of v in element (x,y)

}

double & StaggeredGrid::v( int i, int j) {
	// Part 3/3
	// Access value of v in element (x,y) 
} 	

int StaggeredGrid::vIBegin () const {
	// First valid index for v in x direction 
}

int StaggeredGrid::vIEnd ()	const {
	// One after last valid index for v in x direction 
}

int StaggeredGrid::vJBegin () const {
	// First valid index for v in y direction 
}

int StaggeredGrid::vJEnd ()	const {
	// One after last valid index for v in y direction 
}