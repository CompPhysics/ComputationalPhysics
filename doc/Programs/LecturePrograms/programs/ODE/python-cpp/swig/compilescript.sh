#!/bin/bash

swig -python -c++ -I.. pendelum.i
g++ -fPIC -O -I.. -I/usr/include/python2.4/ -c ../pendelum.cpp pendelum_wrap.cxx
g++ -shared -o _pendelum.so pendelum.o pendelum_wrap.o
