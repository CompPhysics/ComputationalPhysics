# coding=utf-8
# Code to test and demonstrate ODE solvers in the library.
# In order to use rk2 or euler instead of rk4, just change rk4 -> rk2
# or rk4 -> euler in the library function calls
#Written by Kyrre Ness Sjøbæk

import computationalLib, numpy
lib = computationalLib.pylib()

def derivsSHM(y,t):
    """
    Example derivs for solving SHM equation
    x'' + x = 0
    which can be expressed as two coupled diff eq's
    for the variables x and v:
    x'  =  v 
    v' = - x

    Input:
     - Array y[] which contains the values for x, v at time t
     - Array t, the time to evaluate at (useful for external forces etc.)

    I have chosen the following convention:
     - y[0] = x(t)
     - y[1] = v(t)
     
    Output:
     - Array containing the derivatives at point t,
       same convention as with y[]:
       return[0] = x'(t)
       return[1] = v'(t)
    """
    return numpy.array([y[1],-y[0]])

y0 = [1.0,0.0]
t0 = 0.0
te = 10.0
N = 100

#Run and save results in A
A = lib.rk4(y0,t0,te,N, derivsSHM)

for i in xrange (N):
    print str(A[0][i]) + "\t" + str(A[1][i])

#Run and write results to file (GNUPLOT it!)
lib.rk4(y0,t0,te,N, derivsSHM, filename="out.dat")
