#Program which solves the 2+1-dimensional wave equation by a finite difference scheme
#and creates a movie of the time evolution with the scitools.easyviz library.
#The Matlab backend gives a prettier movie than the default gnuplot backend.
#It seems, however that the show=False option does not work for Matlab, hence
#generating the movie takes much more time.
#Written by Magnar K. Bugge

from numpy import *
from scitools.easyviz import *

#Define the grid
N = 31
h = 1.0 / (N-1)
dt = .0005
t_steps = 10000
x,y = ndgrid(linspace(0,1,N),linspace(0,1,N),sparse=False)

alpha = dt**2 / h**2

#Initial conditions with du/dt = 0
u = sin(x*pi)*cos(y*pi-pi/2)
u_old = zeros(u.shape,type(u[0,0]))
for i in xrange(1,N-1):
    for j in xrange(1,N-1):
        u_old[i,j] = u[i,j] + (alpha/2)*(u[i+1,j] - 4*u[i,j] + u[i-1,j] + u[i,j+1] + u[i,j-1])
u_new = zeros(u.shape,type(u[0,0]))

#We don't necessarily want to plot every time step. We plot every n'th step where
n = 100
plotnr = 0

#Iteration over time steps
for k in xrange(t_steps):
    for i in xrange(1,N-1): #1 - N-2 because we don't want to change the boundaries
        for j in xrange(1,N-1):
            u_new[i,j] = 2*u[i,j] - u_old[i,j] + alpha*(u[i+1,j] - 4*u[i,j] + u[i-1,j] + u[i,j+1] + u[i,j-1])

    #Plot
    if k % n == 0:
        plotnr += 1
        mesh(x,y,u_new,hardcopy='frame%04d.png'%plotnr,show=False,axis=[0,1,0,1,-1,1])

    #Prepare for next time step by manipulating pointers
    temp = u_new
    u_new = u_old
    u_old = u
    u = temp

#Make movie
movie('frame*.png',encoder='convert', output_file='movie.gif', fps=10)
