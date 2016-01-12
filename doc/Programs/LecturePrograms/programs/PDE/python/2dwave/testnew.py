#!/usr/bin/env python
# This script reads in data from file with the solutions of the
# 2dim wave function. The data are organized as 
# time 
# l, i, j, u(i,j)   where k is the time index t_l, i refers to x_i and j to y_j
# At the end it converts a series of png files to a movie
# file movie.gif.  You can run this movie file using the imagemagick
# software animate as   - animate movie.gif et voila', Hollywood next

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
ifile = open(inputfilename, 'r')

uarray = {}
for line in ifile:
    elements = line.split()
    if len(elements) < 4: continue
    l, i, j, uvalue = elements
    l = int(l); i = int(i)
    j = int(j); uvalue = float(uvalue)
    if not l in uarray:
        uarray[l] = {}
    if not i in uarray[l]:
        uarray[l][i] = {}
    if not j in uarray[l][i]:
        uarray[l][i][j] = {}
    uarray[l][i][j] = uvalue
ifile.close()

# Fixed Lengths  used in other function to set up the grids.
Lx = 1; nx = len(uarray[1]);
Ly = 1; ny = len(uarray[1][0]);
ntime = len(uarray);
print nx, ny, ntime
x,y = ndgrid(linspace(0,Lx,nx),linspace(0,Ly,ny),sparse=False)

plotnr = 0

u = zeros([nx, ny])
# Loop over time steps
for l in xrange(1, ntime):
    for i in xrange(0,nx):
            for j in xrange(0,ny):
                u[i,j] = uarray[l][i][j]
    plotnr += 1
    mesh(x,y,u,hardcopy='frame%04d.png'%plotnr,show=False,axis=[0,1,0,1,-1,1])

#Make movie
movie('frame*.png',encoder='convert', output_file='movie.gif', fps=10)
cmd = 'animate movie.gif'
os.system(cmd)






