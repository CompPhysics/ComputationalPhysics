from numpy import *
import sys, os

try:
    outfilename = sys.argv[1]
except:
    print "Usage of this script", sys.argv[0], "outfile", sys.argv[1]; sys.exit(1)

# Read file with data
ofile = open(outfilename, 'w')
N = 30
t_steps = 10000
n = 100
l = 0
ofile.write('%d %d\n' % (N,n))
#Iteration over time steps
for k in xrange(t_steps):
    #Print out
    if k % n == 0:
       l += 1
       ofile.write('%g\n' % (k))
       for i in xrange(0,N):
           x = float(i)/N
           for j in xrange(0,N):
               y = float(j)/N
               u = cos(k*pi*sqrt(5))*sin(x*pi)*cos(y*pi*2.0)
               ofile.write('%d %d %d %12.5e\n' % (l,i,j,u))

