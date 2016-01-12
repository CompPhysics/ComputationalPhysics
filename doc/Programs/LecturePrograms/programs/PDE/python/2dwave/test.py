#Program which solves the 2+1-dimensional wave equation by a finite difference scheme
#and creates a movie of the time evolution with the scitools.easyviz library.

from numpy import *
from scitools.easyviz import *

#Define the grid
N = 31
h = 1.0 / (N-1)
dt = .0005
t_steps = 10000
x,y = ndgrid(linspace(0,1,N),linspace(0,1,N),sparse=False)

n = 100
plotnr = 0

#Iteration over time steps
for k in xrange(t_steps):
    u = cos(k*pi*sqrt(5))*sin(x*pi)*cos(y*pi*2.0)
    #Plot
    if k % n == 0:
        plotnr += 1
        mesh(x,y,u,hardcopy='frame%04d.png'%plotnr,show=False,axis=[0,1,0,1,-1,1])

#Make movie
movie('frame*.png',encoder='convert', output_file='movie.gif', fps=10)
