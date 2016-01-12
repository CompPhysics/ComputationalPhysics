#!/usr/bin/env python
# This script reads in data from file with the solutions of the
# 2dim wave function. The data are organized as 
# time 
# l, i, j, u(i,j)   where k is the time index t_l, i refers to x_i and j to y_j
# At the end it converts a series of png files to a movie
# file movie.gif. You can run this movie file using the ImageMagick
# software animate as - animate movie.gif et voila', Hollywood next

# It creates a movie of the time evolution with the scitools.easyviz library.
# To fetch this addition to python go to the link 
# http://code.google.com/p/scitools/wiki/Installation
# This additional tool is the same as that used in INF1100 and should
# be installed on most machines.

from numpy import *
from scitools.easyviz import *
import sys, os

try:
    inputfilename = sys.argv[1]
except:
    print "Usage of this script", sys.argv[0], "inputfile"; sys.exit(1)

# Read file with data
ifile = open(inputfilename)

# Fixed Lengths used in other function to set up the grids.
Lx = 1; nx = 30;
Ly = 1; ny = 30;
ntime = 100;

x, y = ndgrid(linspace(0, Lx, nx), linspace(0, Ly, ny), sparse=False)

plotnr = 0
u = zeros([nx, ny])
# Loop over time steps
for l_ind in xrange(1, ntime + 1):
    for i_ind  in range(0, nx):
        for j_ind in range(0, ny):
            elements = []
            while len(elements) < 4:
                elements = ifile.readline().split()
            l, i, j, value = elements
            if l_ind != int(l): 
                raise IndexError, 'l_ind=%d, l=%d -> l_ind != l' %(l_ind, int(l))
            u[int(i), int(j)] = float(value)
    plotnr += 1
    surf(x, y, u, hardcopy='frame%04d.png' %plotnr, show=True,
         axis=[0, 1, 0, 1,- 1, 1], legend='timestep=%03d' %l_ind,
         shading='interp')
# surf er bedre enn mesh?
# dropper akser ved å sette shading='interp'
# mencoder gir avi-fil
# vbitrate setter opp bitraten

# Make movie
movie('frame*.png', encoder='mencoder', output_file='movie.avi',
      fps=10, vbitrate=2000, qscale=1)
cmd = 'gnome-open movie.avi'
os.system(cmd)


