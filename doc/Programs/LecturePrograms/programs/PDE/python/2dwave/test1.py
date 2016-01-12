#Program which solves the 2+1-dimensional wave equation by a finite difference scheme
#and creates a movie of the time evolution with the scitools.easyviz library.

from numpy import *
from scitools.easyviz import *
import os
from Numeric import *

#Define the grid
t_steps = 10000
Lx = 1; nx = 30; Ly = 1; ny = 30
dx = Lx/float(nx); x =  zeros(nx+1,Float)
dy = Ly/float(ny); y =  zeros(ny+1,Float)
for i in xrange(0,nx):
    x[i] = i*dx
for j in xrange(0,ny):
    y[j] = i*dy
#u = zeros((nx+1,ny+1),Float)

n = 100
plotnr = 0


# Loop over time steps
for k in xrange(t_steps):
    u = cos(k*pi*sqrt(5))*sin(x*pi)*cos(y*pi*2.0)
    # Plot
    if k % n == 0:
        plotnr += 1
        mesh(x,y,u,hardcopy='frame%04d.png'%plotnr,show=False,axis=[0,1,0,1,-1,1])

#Make movie
movie('frame*.png',encoder='convert', output_file='movie.gif', fps=10)
cmd = 'animate movie.gif'
os.system(cmd)
