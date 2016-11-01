For linux/ubuntu users, you need to install two packages (alternatively use the synaptic package manager)
o  sudo apt-get install libopenmpi-dev
o  sudo apt-get install openmpi-bin

For OS X users, install brew (after having installed xcode and gcc, needed for the 
gfortran compiler of openmpi) and then run
o  brew install open-mpi

When running an executable (code.x), run as
o mpirun -n 10 ./code.x
where -n indicates the number of processes, 10 here.

In order to compile, compile as for c++
o mpic++ -O3 -o <executable>  <code.cpp>
and as
o mpif90 -O3 -o <executable>  <code.f90>

With openmpi installed, when using Qt, add to your .pro file the instructions at
http://dragly.org/2012/03/14/developing-mpi-applications-in-qt-creator/

You may need to tell Qt where opempi is stored.

For the machines at the computer lab, openmpi is located  at /usr/lib64/openmpi/bin
Add to your .bashrc file the following
export PATH=/usr/lib64/openmpi/bin:$PATH 


There are several files here
o  MPIsing uses the Ran2 RNG from Numerical Recipes and plain c++ array declaration
o  Paraising.cpp uses the Mersenne twister generator and armadillo to allocate arrays
o  MPIisingOwnclass.cpp	use the Mersenne RNG and my own vectormatrixclass.cpp
o  MPIising.f90 is the corresponding fortran version that uses the xorshift1024* RNG. It needs and input file input.dat
o  MPIintegration.cpp   performs a simple integration in parallel using the trapezoidal rule 
