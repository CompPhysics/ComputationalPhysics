#!/usr/bin/env python
# coding=utf-8
#General RK4 implementation
#Written by Kyrre Ness Sjøbæk

import numpy

def rk4_step(y,t,h,deriv):
    """
    General RK4 stepper for
    N coupled differential eq's

    Input:
     - y:      Array containing the y(t)
     - t:      Which time are we talking about?
     - h:      Stepsize
     - deriv:  Function that returns an array
               containing dy/dt for each y at time t,
               and takes as arguments an y-array, and time t.
    """
    
    k1 = h*deriv(y,t);
    k2 = h*deriv(y+k1/2.0,t+h/2.0)
    k3 = h*deriv(y+k2/2.0,t+h/2.0)
    k4 = h*deriv(y+k3,t+h)

    return y + (k1 + 2*(k2+k3) + k4)/6.0

def rk4(y0, t0, te, N, deriv, filename=None):
    """
    General RK4 driver for
    N coupled differential eq's,
    fixed stepsize

    Input:
     - y0:       Vector containing initial values for y
     - t0:       Initial time
     - te:       Ending time
     - N:        Number of steps
     - deriv:    See rk4_step
     - filename: Optional, use if you want to write
                 data to file at each step.
                 Format used:
                 t y[0] y[1] ... (%10.15E)

    Output:
    If filename=None, return tuple containing:
     - time:  Array of times at which it has iterated over
     - yout:  N*len(y0) numpy array containing y for each timestep
    If filename specified, None is returned.
    
     """

    h = (te-t0)/float(N)
    t = t0;

    if filename == None:
        #Setup arrays
        time = numpy.zeros(N);
        yout = []
        #Inital values
        yout.append(y0);
        time[0] = t0;
        t = t0;

        #Loop over timesteps
        for i in xrange(1,N):
            yout.append(rk4_step(yout[i-1],t,h,deriv));
            t = t0 + h*i;
            time[i] = t;

        return (time,yout)
    else:
        ofile = open(filename,'w')
        #Format string used for output file
        ostring = "%20.8E " + ("%20.8E "*len(y0)) + "\n"

        #Initial values
        y = y0
        t = t0

        foo = [t]; foo[1:] = y;
        ofile.write(ostring % tuple(foo))
        
        while (t < te):
            y = rk4_step(y,t,h,deriv)
            t +=h
            
            foo = [t]; foo[1:] = y;
            ofile.write(ostring % tuple(foo))

        ofile.close()
        return None

# *** EXAMPLE DERIVS: ***

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

def derivThrow(y,t):
    """
    Example derivs for solving a 2D-throwing differential equation.
    In this example gravity is given by the Newtonian expression
    and there is an air resistance that depends linearly on the
    height of the object until it reaches zero.
    
    The y-direction is upward, so the to coupled second order equations are:
    	x''(t) = -(c-a*y(t))*(x'(t))^2
	y''(t) = -g/(r + y(t))^2 -(c - a*y(t))(y'(t))^2
    The constants used here are in calculated in units so that the
    distances are given in kilometers (above the earth for y),
    the time given in hours and the velocities in km/h.
    But this can easily be changed by changing the values of the constants.
    
    The convention for the y-vector is:
	y[0] = x(t)
	y[1] = x'(t)
	y[2] = y(t)
	y[3] = y'(t)
    The values of the constants a and c in the expression
    for the drag coefficient are rather randomly chosen for a rocketlike
    low-resistance object. a was chosen so that the resistance at 1000 km would be zero.
    The if-branch for the resistance was added to ensure that the air resistance
    would never be a force in the direction of the motion, and to ensure sensible
    values slightly below the earth (an obvious improvement to this program would be
    to end the simulation when the object came back down to the earths surface.)
    
    Usage example:
    y0 = numpy.array([0,200,0,300])
    t0 = 0
    te = 2
    N  = 10000
    rk4(y0, t0, te, N, derivThrow, filename='resultat.txt')

    This example was contributed by Marit Sandstad
    """
    a = 1.9E-4
    g = (6.67E-17 * 5.9736E24)
    c = 0.19
    r = 6371 
    b = c - a * y[2]
    if b < 0:
	b = 0
    elif b > 0.19:
	b = 0.19
    print b
    deriv = numpy.zeros(4)
    deriv[0] = y[1]
    deriv[1] = - b * (y[1]**2)
    deriv[2] = y[3]
    deriv[3] = - g / ((y[2] + r)**2) - b * (y[3]**2)

    return deriv
