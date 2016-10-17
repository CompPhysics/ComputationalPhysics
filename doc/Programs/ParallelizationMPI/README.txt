For linux/ubuntu users, you need to install two packages (alternatively use the synaptic package manager)
o  sudo apt-get install libopenmpi-dev
o  sudo apt-get install openmpi-bin

For OS X users, install brew (after having installed xcode and gcc, needed for the 
gfortran compiler of openmpi) and then run
o  brew install open-mpi

When running an executable (code.x), run as
o mpirun -n 10 ./code.x

where -n indicates the number of processes, 10 here.

With openmpi installed, when using Qt, add to your .pro file the instructions at
http://dragly.org/2012/03/14/developing-mpi-applications-in-qt-creator/

You may need to tell Qt where opempi is stored.

For the machines at the computer lab, openmpi is located  at /usr/lib64/openmpi/bin
Add to your .bashrc file the following
export PATH=/usr/lib64/openmpi/bin:$PATH 

