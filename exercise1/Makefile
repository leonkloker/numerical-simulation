CXX=g++         # C++ compiler
CXXFLAGS=-Wall  # compiler flags

# source files
SRC=main.cpp computation.cpp discretization/discretization.cpp discretization/central_differences.cpp discretization/staggered_grid.cpp
   discretization/donor_cell.cpp output_writer/output_writer.cpp output_writer/output_writer_text.cpp output_writer/output_writer_paraview.cpp
   settings/settings.cpp solver/gauss_seidel.cpp solver/sor.cpp solver/pressure_solver.cpp storage/array2d.cpp storage/fieldvariable.cpp


all: debug release

debug:
	$(CXX) $(SRC) $(CXXFLAGS) -g -o test2

release:
	$(CXX) $(SRC) $(CXXFLAGS) -Ofast -DNDEBUG -o test2_release

clean:
	rm test2 test2_release
