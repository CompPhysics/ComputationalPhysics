
#!/usr/bin/env python

# This script reads in data from file with the solutions of the
# 2dim wave function. The data are organized as 
# time 
# x1, y1, u(x1,y1)   
# x2, y1, u(x2,y1)  etc
# till the next time step.  An example is given by the
# wavedata.dat file under project 4.
# It uses gnuplot's splot function. 
# At the end it converts a series of png files to a movie
# file movie.gif.  You can run this movie file using the imagemagick
# software animate as   - animate movie.gif et voila', Hollywood next
# At the end, it prepares the body of a latex file and an html report
# as well.

import sys, os
from Numeric import *
import Gnuplot

g = Gnuplot.Gnuplot(persist=1)

# The script is run as python script.py filename.dat

try:
    infilename = sys.argv[1]
except:
    print "Usage of this script", sys.argv[0], "infile", sys.argv[1]; sys.exit(1)

# Read file with data
ifile = open(infilename, 'r')
lines = ifile.readlines()
# Fill in x, y and u(x,y)
x = [] ;  y = []   #  This are empty lists
for line in lines:
    xvalue, yvalue = line.split()
    x.append(float(xvalue)); y.append(float(yvalue))
ifile.close()

# convert to a form that the gnuplot interface can deal with
#d = Gnuplot.Data(x, y, title='data from output file', with='lp')
d = Gnuplot.Data(x, y, u, title='data from output file')
g.xlabel('x')   #  make x label
g.ylabel('y') 
g.zlabel('u(x,y)') 
g.splot(d)                         # plot the data
g.hardcopy(filename="relerror.ps",terminal="postscript", enhanced=1, color=1)

# Now prepare latex file, r in front avoids backslashes being treated
# as control chars in strings
preamb = r"""\documentclass[prc,aps,twocolumn,floatfix]{revtex4}
\usepackage[dvips]{graphicx}
\usepackage{epsfig}
\usepackage{pst-plot}
\usepackage{bm}

   """

title = "Begin report 4"

figure = r"""\begin{figure}[hbtp]
\includegraphics[angle=270, scale=0.35]{relerror.ps}
\caption{Simple plot} 
\label{fig:simpleplot}
\end{figure}

"""


# Dump to file:
filename = 'report_project'
f = file(filename + '.tex', "w")
f.write(preamb)
f.write(r"""\begin{document}""")
f.write(title)
f.write(figure)
f.write("""\end{document}""")
f.close()


#Define the grid
N = 31
h = 1.0 / (N-1)
dt = .0005
t_steps = 10000
x,y = ndgrid(linspace(0,1,N),linspace(0,1,N),sparse=False)



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


#Make movie
movie('frame*.png',encoder='convert', output_file='movie.gif', fps=10)



