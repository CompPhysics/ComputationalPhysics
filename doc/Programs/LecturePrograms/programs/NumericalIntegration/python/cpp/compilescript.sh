#!/bin/bash

swig -python -c++ -I. lib.i
g++ -fPIC -O -I.. -I/usr/include/python2.4/ -I/mn/felt/u2/kyrrens/local/include/ -c lib.cpp lib_wrap.cxx
g++ -shared -o _pylib_cpp.so lib.o lib_wrap.o
