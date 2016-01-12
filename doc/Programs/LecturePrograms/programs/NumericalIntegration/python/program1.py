# -*- coding: utf-8 -*-
#Example of numerical integration with Gauss-Legendre quadrature
#Translated to Python by Kyrre Ness Sjøbæk

import sys
import numpy
from computationalLib import pylib

#Read input
if len(sys.argv) == 1:
    print "Number of integration points:"
    n = int(sys.stdin.readline())
    print "Integration limits (a, b)"
    (a,b) = sys.stdin.readline().split(",")
    a = int(a)
    b = int(b)
elif len(sys.argv) == 4:
    n = int(sys.argv[1])
    a = int(sys.argv[2])
    b = int(sys.argv[3])
else:
    print "Usage: python", sys.argv[0], "N a b"
    sys.exit(0)

#Definitions
m = pylib(inputcheck=False,cpp=False)

def integrand(x):
    return 4./(1. + x*x)

#Integrate with Gaus-legendre!
(x,w) = m.gausLegendre(a,b,n)
int_gaus = numpy.sum(w*integrand(x))

#Final output
print "integration of f(x) 4/(1+x**2) from",a,"to",b,"with",n,"meshpoints:"
print "Gaus-legendre:", int_gaus
print "Trapezoidal:  ", m.trapezoidal(a,b,n, integrand)
print "Simpson:      ", m.simpson(a,b,n,integrand)
