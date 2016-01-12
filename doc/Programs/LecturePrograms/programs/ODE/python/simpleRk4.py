# coding=utf-8
#
# This program solves Newtons equation for a block sliding on
# an horizontal frictionless surface.
# The block is tied to the wall with a spring, so N's eq takes the form:
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
# Note that this is a highly simplifyed rk4 code, intended
# for conceptual understanding and experimentation.
# A more generalized code is also aviable.
#
#Written by Kyrre Ness Sjøbæk

import sys
import numpy, math

#Global variables
ofile = None;
E0    = 0.0

def sim(x_0, v_0, N):
    ts = 0.0
    te = 4*math.pi
    h = (te-ts)/float(N)

    t = ts;
    x = x_0
    v = v_0
    while (t < te):
        kv1 = -h*x
        kx1 = h*v

        kv2 = -h*(x+kx1/2)
        kx2 =  h*(v+kv1/2)

        kv3 = -h*(x+kx2/2)
        kx3 =  h*(v+kv2/2)

        kv4 = -h*(x+kx3/2)
        kx4 =  h*(v+kv3/2)

        #Write the old values to file
        output(t,x,v)

        #Update
        x = x + (kx1 + 2*(kx2+kx3) + kx4)/6
        v = v + (kv1 + 2*(kv2+kv3) + kv4)/6
        t = t+h
        
def output(t,x,v):
    de = 0.5*x**2+0.5*v**2 - E0;
    ofile.write("%15.8E %15.8E %15.8E %15.8E %15.8E\n"\
                %(t, x, v, math.cos(t),de));


#MAIN PROGRAM:

#Get input
if len(sys.argv) == 5:
    ofilename =       sys.argv[1];
    x_0       = float(sys.argv[2])
    v_0       = float(sys.argv[3])
    N         =   int(sys.argv[4])
else:
    print "Usage:", sys.argv[0], "ofilename x0 v0 N"
    sys.exit(0)

#Setup
ofile = open(ofilename, 'w')
E0    = 0.5*x_0**2+0.5*v_0**2

#Run simulation
sim(x_0,v_0,N)

#Cleanup
ofile.close() 
