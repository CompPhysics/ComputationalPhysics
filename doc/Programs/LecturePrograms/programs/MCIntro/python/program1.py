# coding=utf-8
#Example of brute force Monte Carlo integration
#Translated to Python by Kyrre Ness Sjøbæk

import math, numpy, sys

def func(x):
    """Integrand"""
    return 4/(1.0+x*x)

#Read in number of samples
if len(sys.argv) == 2:
    N = int(sys.argv[1])
else:
    print "Usage:",sys.argv[0],"N"
    sys.exit(0)

#Evaluate the integral using a crude MonteCarlo
mc_sum  = 0.0
mc_sum2 = 0.0
for i in xrange(N):
    fx=func(numpy.random.random())
    mc_sum  += fx
    mc_sum2 += fx*fx

#Fix the statistics
mc_sum  /= float(N)
mc_sum2 /= float(N)

variance = mc_sum2 - mc_sum**2

print "Integral=",mc_sum,", variance=",variance
    
