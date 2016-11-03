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

You can also use mpicxx instead of mpic++ and mpiexec instead of mpirun.

With openmpi installed, when using Qt, add to your .pro file the instructions at
http://dragly.org/2012/03/14/developing-mpi-applications-in-qt-creator/

You may need to tell Qt where opempi is stored. This can be done by clicking on the project link
and scroll down to build environment. Here you update the PATH. In order to find where mpic++ is use
o which mpcic++   or which mpicxx
It can for example be found in /usr/local/bin. Add this to the PATH. 

For the machines at the computer lab, openmpi is located  at /usr/lib64/openmpi/bin
Add to your .bashrc file the following
export PATH=/usr/lib64/openmpi/bin:$PATH 


There are several files here
o  MPIsing uses the Ran2 RNG from Numerical Recipes and plain c++ array declaration
o  Paraising.cpp uses the Mersenne twister generator and armadillo to allocate arrays
o  MPIisingOwnclass.cpp	use the Mersenne RNG and my own vectormatrixclass.cpp to allocate arrays
o  MPIising.f90 is the corresponding fortran version that uses the xorshift1024* RNG. It needs and input file input.dat
o  MPIintegration.cpp   performs a simple integration in parallel using the trapezoidal rule 
o  MPIHelloworld.cpp is pretty obvious!

For running on SMAUG (our local cluster), go to http://comp-phys.net/ and click on the link internals and click on computing cluster. To get access to Smaug, you will need to send us an e-mail with your name, UiO username, phone number, room number and affiliation to the research group. State that in the subjectline that your are a participants of the course FYS3150/4150 and write that you are aware that your account has to be deleted when the semester is over. In return, you will receive a password you may use to access the cluster.

Here follows a simple recipe

   log in as ssh -l username tid.uio.no
   ssh username@fyslab-compphys
In the folder

    shared/guides/starting_jobs 
you will find a simple example on how to set up a job and compile and run. This files are write protected. Copy them to your own folder and compile and run there. For more information see the readme file under the program folder.
