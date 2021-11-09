#include "staggered_grid.h"
 
StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) : nCells_(nCells), meshWidth_(meshWidth), 
u_(FieldVariable({nCells_[0]+2, nCells_[1]+2}, {0, -0.5*meshWidth[1]}, meshWidth_)), 
v_(FieldVariable({nCells_[0]+2, nCells_[1]+2}, {-0.5*meshWidth[0],0}, meshWidth_)),
p_(FieldVariable({nCells_[0]+2, nCells_[1]+2}, {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth_)),
rhs_(FieldVariable({nCells_[0]+2, nCells_[1]+2}, {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth_)),
f_(FieldVariable({nCells_[0]+2, nCells_[1]+2}, {0, -0.5*meshWidth[1]}, meshWidth_)),
g_(FieldVariable({nCells_[0]+2, nCells_[1]+2}, {-0.5*meshWidth[0],0}, meshWidth_)){}

double StaggeredGrid::dx () const {
	return meshWidth_[0];
}

double StaggeredGrid::dy () const {
	return meshWidth_[1];
}

double & StaggeredGrid::f(int i, int j){
	return f_(i,j);
} 		

double & StaggeredGrid::g(int i, int j){
	return g_(i,j);
} 		

const std::array< double, 2 > StaggeredGrid::meshWidth() const {
	return meshWidth_;
}

const std::array< int, 2 > StaggeredGrid::nCells() const {
	return nCells_;
}

const FieldVariable & StaggeredGrid::p() const {
	const FieldVariable& p_ref = p_;
	return p_ref; 
}
double StaggeredGrid::p(int i,int j) const {
	return p_(i,j);
}

double & StaggeredGrid::p(int i, int j) {
	return p_(i,j);
}

int StaggeredGrid::pIBegin() const {
	return 0;
}

int StaggeredGrid::pIEnd () const {
	return nCells_[0]+2;
}

int StaggeredGrid::pJBegin() const {
	return 0;
}

int StaggeredGrid::pJEnd() const {
	return nCells_[1]+2;
}

double & StaggeredGrid::rhs(int	i, int j) {
	return rhs_(i,j);
}

const FieldVariable & StaggeredGrid::u () const{
	const FieldVariable& u_ref = u_;
	return u_ref; 
}

double StaggeredGrid::u (int i,	int j) const {
	return u_(i,j);
}

double & StaggeredGrid::u (int i, int j) {
	return u_(i,j);
} 		

int StaggeredGrid::uIBegin () const {
	return 0;
}

int StaggeredGrid::uIEnd () const {
	return nCells_[0]+1;
}

int StaggeredGrid::uJBegin () const {
	return 0;
}

int StaggeredGrid::uJEnd () const {
	return nCells_[1]+2;
}

const FieldVariable & StaggeredGrid::v () const {
	const FieldVariable& v_ref = v_; 
	return v_ref;
}

double StaggeredGrid::v (int i,	int j) const {
	return v_(i,j);
}

double & StaggeredGrid::v(int i, int j) {
	return v_(i,j); 
} 	

int StaggeredGrid::vIBegin () const {
	return 0;
}

int StaggeredGrid::vIEnd () const {
	return nCells_[0]+2;;
}

int StaggeredGrid::vJBegin () const {
	return 0;
}

int StaggeredGrid::vJEnd () const {
	return nCells_[1]+1;
}
