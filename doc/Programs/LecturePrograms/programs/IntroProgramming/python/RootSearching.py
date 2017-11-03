"""
Bisection, Secant & Newton Raphson Method.
"""

import math
"""
* Variable Description:
*
* f			:	Given function
* f_		:	Derivative of f
* [a, b]	:	End point values
* TOL		:	Tolerance
* NMAX		:	Max number of iterations
"""


def bisection(f, a, b, TOL=0.001, NMAX=100):
	"""
	Takes a function f, start values [a,b], tolerance value(optional) TOL and
	max number of iterations(optional) NMAX and returns the root of the equation
	using the bisection method.
	"""
	n=1
	while n<=NMAX:
		c = (a+b)/2.0
		print "a=%s\tb=%s\tc=%s\tf(c)=%s"%(a,b,c,f(c))
		if f(c)==0 or (b-a)/2.0 < TOL:
			return c
		else:
			n = n+1
			if f(c)*f(a) > 0:
				a=c
			else:
				b=c
	return False

def secant(f,x0,x1, TOL=0.001, NMAX=100):
	"""
	Takes a function f, start values [x0,x1], tolerance value(optional) TOL and
	max number of iterations(optional) NMAX and returns the root of the equation
	using the secant method.
	"""
	n=1
	while n<=NMAX:
		x2 = x1 - f(x1)*((x1-x0)/(f(x1)-f(x0)))
		if x2-x1 < TOL:
			return x2
		else:
			x0 = x1
			x1 = x2
	return False

def newtonraphson(f, f_, x0, TOL=0.001, NMAX=100):
	"""
	Takes a function f, its derivative f_, initial value x0, tolerance value(optional) TOL and
	max number of iterations(optional) NMAX and returns the root of the equation
	using the newton-raphson method.
	"""
	n=1
	while n<=NMAX:
		x1 = x0 - (f(x0)/f_(x0))
		if x1 - x0 < TOL:
			return x1
		else:
			x0 = x1
	return False

if __name__ == "__main__":
	
	def func(x):
	"""
	Function x^3 - x -2
	We will calculate the root of this function using different methods.
	"""
		return math.pow(x,3) - x -2

	def func_(x):
	"""
	Derivate of the function f(x) = x^3 - x -2
	This will be used in Newton-Rahson method.
	"""
		return 3*math.pow(x,2)-1
	
	#Invoking Bisection Method
	res = bisection(func,1,2).
	print res
	
	#Invoking Secant Method
	res = bisection(func,1,2).
	print res
	
	#Invoking Newton Raphson Method
	res = newtonraphson(func,func_,1)
	print res
