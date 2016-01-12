#!/bin/bash
swig -python -c++ lib.i
root=`python -c 'import sys; print sys.prefix'`
ver=`python -c 'import sys; print sys.version[:3]'`
#g++ -O -I. -I$root/include/python$ver -c lib.cpp lib_wrap.cxx
g++ -fPIC -I. -I/local/include/python2.4 -c lib.cpp lib_wrap.cxx
g++ -shared -o _lib.so lib.o lib_wrap.o