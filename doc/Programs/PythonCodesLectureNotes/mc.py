from  matplotlib import pyplot as plt
from math import exp, acos, log10
import numpy as np
from sympy import Symbol, integrate, exp, oo
import random


# function for the trapezoidal rule
def TrapezoidalRule(a,b,f,n):
   h = (b-a)/float(n)
   s = 0
   x = a
   for i in range(1,n,1):
       x = x+h
       s = s+ f(x)
   s = 0.5*(f(a)+f(b))+s
   return h*s
# function to perform the Monte Carlo calculations
def MonteCarloIntegration(f,n):
    sum = 0
# Define the seed for the rng
    random.seed()    
    for i in range (1, n, 1):
        x = random.random()
        sum = sum +f(x)
    return sum/n

#  function to compute
def function(x):
    return 4/(1+x*x)

# Integration limits for the Trapezoidal rule
a = 0.0; b = 1.0
# define x as a symbol to be used by sympy
x = Symbol('x')
# find result from sympy
#exact = integrate(function(x), (x, a, b))
exact = acos(-1.0)
# set up the arrays for plotting the relative error
log10n = np.zeros(6); Trapez = np.zeros(6); MCint = np.zeros(6);
# find the relative error as function of integration points
for i in range(1, 6):
    npts = 10**(i+1)
    log10n[i] = log10(npts)
    Trapez[i] = log10(abs((TrapezoidalRule(a,b,function,npts)-exact)/exact))
    MCint[i] = log10(abs((MonteCarloIntegration(function,npts)-exact)/exact))
plt.plot(log10n, Trapez ,'b-',log10n, MCint,'g-')
plt.axis([1,6,-14.0, 0.0])
plt.xlabel('$\log_{10}(n)$')
plt.ylabel('Relative error')
plt.title('Relative errors for Monte Carlo integration and Trapezoidal rule')
plt.legend(['Trapezoidal rule', 'Brute force Monte Carlo integration'], loc='best') 
plt.savefig('mcintegration.pdf')
plt.show()
