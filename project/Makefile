CXX=g++         # C++ compiler
CXXFLAGS=-Wall  # compiler flags

# source files
SRC=src/main.cpp src/computation_parallel.cpp src/discretization/discretization.cpp src/discretization/central_differences.cpp src/discretization/staggered_grid.cpp src/discretization/donor_cell.cpp src/output_writer/output_writer.cpp src/output_writer/output_writer_text.cpp src/output_writer/output_writer_paraview_parallel.cpp src/output_writer/output_writer_text_parallel.cpp src/settings/settings.cpp src/solver/sor_parallel.cpp src/solver/cg_parallel.cpp src/solver/pressure_solver_parallel.cpp src/temperature/temperature.cpp src/storage/array2d.cpp src/storage/fieldvariable.cpp src/partitioning/partitioning.cpp


all: debug release

debug:
	$(CXX) $(SRC) $(CXXFLAGS) -g -o project

release:
	$(CXX) $(SRC) $(CXXFLAGS) -Ofast -DNDEBUG -o project

clean:
	rm test2 project
