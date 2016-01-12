# coding=utf-8
# Code to test and demonstrate the adaptive-step rk4 found in the library.
#Written by Kyrre Ness Sjøbæk

import computationalLib, numpy
lib = computationalLib.pylib()

import math

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

#only return lowh once
lowhChanged = False; 

def errorfuncSHM_energy (yOld, yNew, yErr, h, t):
    """
    Example errorfunc for SHM, using conservation of energy
    as the criterion.

    Doubles or halves the stepsize.
    """
    tolerance = 0.001      # Accept if the error is below 0.1%
    lowtol    = 0.0005     # If error below 0.05%, double stepsize
    lowh      = 0.0001 # Lower bound for h - important to implement so
                          # we don't end up in an infinite loop

    energyNow = 0.5*(yNew[0]**2 + yNew[1]**2)
    
    test      = abs(energyNow-0.5)/0.5; # I know (analyticaly, initial conds.)
                                        # that E = 0.5

    #print yNew, t, h, test

    global lowhChanged

    if test < tolerance:
        #Accept. Are the error "to low"?
        if test < lowtol:
            if lowhChanged:
                print "Restored at time %.7E" % t
                lowhChanged = False
            return 2*h;
        else:
            return 0.0
    else:
        #To big error, reduce stepsize
        hnew = h/2.0
        if hnew < lowh:
            print "**** h = lowh at time %.7E - solution may be innacurate! ****" % t
            if lowhChanged:
                return 0.0
            else:
                lowhChanged = True
                return lowh;
        else:
            return hnew
        
    


y0 = [1.0,0.0]
t0 = 0.0
te = 50.0
h0 = 0.1

#Run and save results in A:

#Use energy conservation as the stepsize requirement
#A = lib.rk4Adaptive(y0,t0,te,h0, derivsSHM, errorfuncSHM_energy)

#Use fixed stepsize
#A = lib.rk4Adaptive(y0,t0,te,h0, derivsSHM, lambda yOld, yNew, yErr, h, t: 0.0)

#Use fixed tolerance:
#lib.yTol = numpy.array([1,1])*1e-6 #epsilon = 10^-6
#A = lib.rk4Adaptive(y0,t0,te,h0, derivsSHM, lib.rk4Adaptive_stepsizeControl2)

#Write results to screen
#for i in xrange (len(A[0])):
#    print str(A[0][i]) + "\t" + str(A[1][i])

#Run and write results to file (GNUPLOT it!):
lib.rk4Adaptive(y0,t0,te,h0, derivsSHM, errorfuncSHM_energy, filename="out.dat")
