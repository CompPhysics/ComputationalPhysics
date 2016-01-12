#!/bin/bash

#Compile swig interface
swig -python -c++ -I. ising2dim_backend.i
root=`python -c 'import sys; print sys.prefix'`
ver=`python -c 'import sys; print sys.version[:3]'`
g++ -fPIC -O -I. -I/usr/include/python2.4/ -c ising2dim_backend.cpp lib.cpp ising2dim_backend_wrap.cxx
g++ -shared -o _ising2dim_backend.so ising2dim_backend.o ising2dim_backend_wrap.o lib.o