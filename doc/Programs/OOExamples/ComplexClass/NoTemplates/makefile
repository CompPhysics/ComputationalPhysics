# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3 -std=c++11
# Here we define the library functions we need (not needed for this particular code though)
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= complex
${PROG} :	   	main.o  mycomplex.o 
			${CPPflags} main.o mycomplex.o  ${LIB} -o ${PROG}

main.o :			main.cpp 
		        	${CPPflags} -c main.cpp

mycomplex.o :			mycomplex.cpp
		        	${CPPflags} -c mycomplex.cpp

