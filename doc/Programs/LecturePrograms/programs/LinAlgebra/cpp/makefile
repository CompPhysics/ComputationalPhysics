GCC= c++ 
PROG= invers

${PROG} :	   ${PROG}.o lib.o
		   ${GCC} ${PROG}.o lib.o -o ${PROG}

${PROG}.o :	   ${PROG}.cpp
		   ${GCC} -c ${PROG}.cpp

lib.o :		   lib.cpp lib.h
		   ${GCC} -c lib.cpp

