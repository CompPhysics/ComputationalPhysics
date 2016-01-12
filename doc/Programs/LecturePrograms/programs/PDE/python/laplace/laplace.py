# coding=utf-8
# Solves the 2d Laplace equation using relaxation method
#Written by Kyrre Ness Sjøbæk

import numpy, math

def relax(A, maxsteps, convergence):
    """
    Relaxes the matrix A until the sum of the absolute differences
    between the previous step and the next step (divided by the number of
    elements in A) is below convergence, or maxsteps is reached.

    Input:
     - A: matrix to relax
     - maxsteps, convergence: Convergence criterions

    Output:
     - A is relaxed when this method returns
    """

    iterations = 0
    diff = convergence +1

    Nx = A.shape[1]
    Ny = A.shape[0]
    
    while iterations < maxsteps and diff > convergence:
        #Loop over all *INNER* points and relax
        Atemp = A.copy()
        diff = 0.0
        
        for y in xrange(1,Ny-1):
            for x in xrange(1,Ny-1):
                A[y,x] = 0.25*(Atemp[y,x+1]+Atemp[y,x-1]+Atemp[y+1,x]+Atemp[y-1,x])
                diff  += math.fabs(A[y,x] - Atemp[y,x])

        diff /=(Nx*Ny)
        iterations += 1
        print "Iteration #", iterations, ", diff =", diff;


def boundary(A,x,y):
    """
    Set up boundary conditions

    Input:
     - A: Matrix to set boundaries on
     - x: Array where x[i] = hx*i, x[last_element] = Lx
     - y: Eqivalent array for y

    Output:
     - A is initialized in-place (when this method returns)
    """

    #Boundaries implemented (condensator with plates at y={0,Lx}, DeltaV = 200):
    # A(x,0)  =  100*sin(2*pi*x/Lx)
    # A(x,Ly) = -100*sin(2*pi*x/Lx)
    # A(0,y)  = 0
    # A(Lx,y) = 0

    Nx = A.shape[1]
    Ny = A.shape[0]
    Lx = x[Nx-1] #They *SHOULD* have same sizes!
    Ly = x[Nx-1]


    A[:,0]    =   100*numpy.sin(math.pi*x/Lx)
    A[:,Nx-1] = - 100*numpy.sin(math.pi*x/Lx)
    A[0,:]    = 0.0
    A[Ny-1,:] = 0.0


#Main program

import sys

if len(sys.argv) == 5:
    Nx      = int(sys.argv[1])
    Ny      = int(sys.argv[2])
    maxiter = int(sys.argv[3])
    ofile   =     sys.argv[4]
else:
    print "Usage:", sys.argv[0], "Nx Ny maxiter outfile"
    print "Lx = Ly = 1, Nx and Ny refers to number of INNER meshpoints"
    sys.exit(0)

x = numpy.linspace(0,1,num=Nx+2) #Also include edges
y = numpy.linspace(0,1,num=Ny+2)
A = numpy.zeros((Nx+2,Ny+2))

boundary(A,x,y)
#Remember: as solution "creeps" in from the edges,
#number of steps MUST AT LEAST be equal to
#number of inner meshpoints/2 (unless you have a better
#estimate for the solution than zeros() )
relax(A,maxiter,0.00001)

numpy.savetxt(ofile, A, fmt="%.10E", delimiter=' ')
