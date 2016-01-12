#
# This program solves genereal 2'nd order diff. eq's using RK4.
# General form of diff eq.:
#
#  dx^2/dt^2 + p(t,x,v)*dx/dt + q(t,x,v)*x = r(t,x,v)
#
# This can generally be decomposed to two coupled 2nd order eq's:
# dx^2/dt^2 = dv/dt = r - p*v - q*x
# dx/dt = v
#
# These equations are implemented by the functions "deriv1", "deriv2" etc.
# (for different equation sets, that is, different r,p,and q)
#
# Implemented example 1 (derivs1): Newtons equation for a block sliding on
# an horizontal frictionless surface. The block is tied to
# the wall with a spring, so N's eq takes the form:
#
#  m d^2x/dt^2 = - kx
#
# In order to make the solution dimless, we set k/m = 1.
# This results in two coupled diff. eq's that may be written as:
#
#  dx/dt = v
#  dv/dt = -x
#
# The user has to specify the initial velocity and position,
# and the number of steps. The time interval is fixed to
# t \in [0, 4\pi) (two periods)
#
# This uses a more advanced and more flexible way of doing a RK4
# simulation, inspired by Numerical Recipes and M.H.Jensen.
# 
# Written by Kyrre N. Sjobak (k.n.sjobak@fys.uio.no)
#

import sys
import numpy, math

#Global variables
ofile = None;
E0    = 0.0

Q = A_roof = omega_0 = omega_roof = g = 0.0
method = 0

def sim1(x_0, v_0, N):
    """
    Driver routine for RK4 integration and file output
    for SHM
    """
    E0 = 0.5*x_0**2+0.5*v_0**2
    te = 4*math.pi
    ts = 0.0
    h = (te-ts)/float(N)

    t = ts;
    x = x_0
    v = v_0
    while (t < te):
        #Write the old values to file
        de = 0.5*x**2+0.5*v**2 - E0;
        ofile.write("%15.8E %15.8E %15.8E %15.8E %15.8E\n"\
                    %(t, x, v, math.cos(t),de));

        #Move one step forward
        (x,v) = rk4(x,v,t,h,derivs1)

        #Update time
        t = t+h

def sim2(x_0, v_0, N, te):
    """
    Driver routine for RK4 integration and file output
    for dimensionless damped pendulum
    """
    delta_t = te/float(N)
    delta_t_roof = omega_0*te/float(N)
    
    t   = 0.0;
    t_h = 0.0;
    x   = x_0;
    v   = v_0;
    while (t < te):
        #Write the old values to file
        ofile.write("%15.8E %15.8E %15.8E\n" %(t, x, v));

        #Move one step forward
        (x,v) = rk4(x,v,t_h,delta_t_roof,derivs2)
        
        #Update time
        t   += delta_t      #file-time
        t_h += delta_t_roof #dimless sim-time
    

def rk4(x, v, t, h, deriv):
    """
    Core Rk4 algo for calculating one step ahead.
    This version is limited to 2nd order eq's of the type
    dx^2/dt^2 + p(t,x,v)*dx/dt + q(t,x,v)*x = r(t,x,v)

    Input:
    - x:      x (t)
    - v:      v (t)
    - t:      Initial time (t_n)
    - h:      Stepsize (t_{n+1} = t_n + h)
    - deriv:  Pointer to a function that calculates
              and returns the derivatives of x and v
    Output:
    (x,v):    Tuple containing the calculated
              values of x and v at time t+h
    """

    (dxdt,dvdt) = deriv(t,x,v)
    kv1 = h*dvdt
    kx1 = h*dxdt

    (dxdt,dvdt) = deriv(t+h/2,x+kx1/2,v+kv1/2)
    kv2 = h*dvdt
    kx2 = h*dxdt
    
    (dxdt,dvdt) = deriv(t+h/2,x+kx2/2,v+kv2/2)
    kv3 = h*dvdt
    kx3 = h*dxdt

    (dxdt,dvdt) = deriv(t+h,x+kx3,v+kv3)
    kv4 = h*dvdt
    kx4 = h*dxdt

    x = x + (kx1 + 2*(kx2+kx3) + kx4)/6
    v = v + (kv1 + 2*(kv2+kv3) + kv4)/6

    return(x,v)
    
def derivs1(t,x,v):
    """
    Derivatives for simple harm. osc. (\"sliding block\")
    """

    dxdt = v
    dvdt = -x
    
    return(dxdt, dvdt)

def derivs2(t,x,v):
    """
    Derivates for pendulum with dampening (dimless).
    Here we use x = theta, v = v_roof
    """

    dxdt = v
    if Q != 0:
        dvdt = A_roof*math.cos(omega_roof*t) - v/Q - math.sin(x)
    else:
        dvdt = A_roof*math.cos(omega_roof*t) - math.sin(x)

    return(dxdt,dvdt)



#Get input
if len(sys.argv) >= 2:
    if len(sys.argv) == 6 and int(sys.argv[1]) == 1:
        method    = 1
        ofilename =       sys.argv[2]
        x_0       = float(sys.argv[3])
        v_0       = float(sys.argv[4])
        N         =   int(sys.argv[5])
    elif len(sys.argv) == 12 and int(sys.argv[1]) == 2:
        method    = 2
        ofilename =       sys.argv[2]
        m         = float(sys.argv[3])
        l         = float(sys.argv[4])
        omega     = float(sys.argv[5])
        A         = float(sys.argv[6])
        viscosity = float(sys.argv[7])
        x_0       = float(sys.argv[8])
        v_0       = float(sys.argv[9])
        te        = float(sys.argv[10]) #Final time as a multiple of pi
        N         =   int(sys.argv[11])
    else:
        print "Usage:", sys.argv[0], "method\n ofilename [methodargs]"
        print "where method = 1 for simple harm. osc, 2 for damped pendulum"
        print
        print "Methodargs for method 1: x_0 v_0 N"
        print "Methodargs for method 2: mass length omega amplitude viscosity theta_0 v_0 te N"
        sys.exit(0)

#Setup
ofile = open(ofilename, 'w')


#Run simulation
if method == 1:
    sim1(x_0,v_0,N)
    
elif method == 2:
    te = te*math.pi
    g  = 9.81
    omega_0 = math.sqrt(g/l)
    if viscosity != 0:
        Q = m*g/(omega_0*viscosity)
    else:
        Q = 0
    A_roof = A/(m*g)
    omega_roof = omega/omega_0

    sim2(x_0,v_0,N,te)
    
    
#Cleanup
ofile.close() 
