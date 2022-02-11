Project submission of group Maksim Helmann, Leon Kloker

We implemented a parallelized conjugate gradient method as an additional solver for the pressure Poisson equation 
as well as a parallel temperature solver with arbitrary Dirichlet and Neumann boundary conditions. We provide three test cases 
which are given in the parameters_case1.txt, parameters_case2.txt and parameters_case3.txt files.

Steps to run the code:

- mkdir build
- cd build
- cmake -DCMAKE_BUILD_TYPE=RELEASE ..
- make install
- mpirun -n #even number of processes ./numsim_parallel ../parameters_case1.txt

