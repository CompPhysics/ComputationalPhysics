# coding=utf-8
#Monte Carlo integration in 6 dimensions
#Translated by Kyrre Ness Sjøbæk

import numpy,math
import sys

def integrand(x):
    """Calculates the integrand
    exp(-a*(x1^2+x2^2+...+x6^2)-b*[(x1-x4)^2+...+(x3-x6)^2])
    from the values in the 6-dimensional array x."""
    a = 1.0
    b = 0.5

    x2 = numpy.sum(x**2)
    xy = (x[0]-x[3])**2 + (x[1]-x[4])**2 + (x[2]-x[5])**2
    
    return numpy.exp(-a*x2-b*xy)

#Main program

#Integration limits: x[i] in (-5, 5)
L      = 5.0
jacobi = (2*L)**6

if len(sys.argv) == 2:
    N = int(sys.argv[1])
else:
    print "Usage: python",sys.argv[0],"number of samples"
    sys.exit(0)

print "Running with N = %d..." % N

#Evaluate the integral
#See Numerical Recipes 7.6
sum  = 0.0
sum2 = 0.0
for i in xrange(N):
    #Generate random coordinates to sample at
    x = numpy.array([-L+2*L*numpy.random.random() for j in xrange(6)])

    fx         = integrand(x)
    sum       += fx
    sum2      += fx**2
#Calculate expt. values for fx and fx^2
sum /=float(N)
sum2/=float(N)

#Result
int_mc  = jacobi*sum;
#Gaussian standard deviation
sigma   = jacobi*math.sqrt((sum2-sum**2)/float(N))

#Output
print "Montecarlo result = %10.8E" % int_mc
print "Sigma             = %10.8E" % sigma
