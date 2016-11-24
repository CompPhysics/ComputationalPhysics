
# Code for solving du/dt = ddu/ddx on a rectangular grid
# of size L x (T*dt),
# with with L = 1, u(x,0) = g(x), u(0,t) = u(L,t) = 0
#

import numpy,math
import sys

def forward_step(alpha,u,uPrev,N):
    """
    Steps forward-euler algo one step ahead.
    Implemented in separate function for code-reuse from crank_nicolson()
    """
    
    for x in xrange(1,N+1): #loop from i=1 to i=N
        u[x] = alpha*uPrev[x-1] + (1.0-2*alpha)*uPrev[x] + alpha*uPrev[x+1]

def forward_euler(alpha,u,N,T):
    """
    Implements forward euler sheme, results saved to
    array u
    """

    #Skip boundary elements
    for t in xrange(1,T):
        forward_step(alpha,u[t],u[t-1],N)

def tridiag(alpha,u,N):
    """
    Tridiagonal gaus-eliminator, specialized to diagonal = 1+2*alpha,
    super- and sub- diagonal = - alpha
    """
    d = numpy.zeros(N) + (1+2*alpha)
    b = numpy.zeros(N-1) - alpha

    #Forward eliminate
    for i in xrange(1,N):
        #Normalize row i (i in u convention):
        b[i-1] /= d[i-1];
        u[i] /= d[i-1] #Note: row i in u = row i-1 in the matrix
        d[i-1] = 1.0
        #Eliminate
        u[i+1] += u[i]*alpha
        d[i] += b[i-1]*alpha
    #Normalize bottom row
    u[N] /= d[N-1]
    d[N-1] = 1.0

    #Backward substitute
    for i in xrange(N,1,-1): #loop from i=N to i=2
        u[i-1] -= u[i]*b[i-2]
        #b[i-2] = 0.0 #This is never read, why bother...
        
def backward_euler(alpha,u,N,T):
    """
    Implements backward euler scheme by gaus-elimination of tridiagonal matrix.
    Results are saved to u.
    """
    for t in xrange(1,T):
        u[t] = u[t-1].copy()
        tridiag(alpha,u[t],N) #Note: Passing a pointer to row t, which is modified in-place

def crank_nicolson(alpha,u,N,T):
    """
    Implents crank-nicolson scheme, reusing code from forward- and backward euler
    """
    for t in xrange(1,T):
        forward_step(alpha/2,u[t],u[t-1],N)
        tridiag(alpha/2,u[t],N)

def g(x):
    """Initial condition u(x,0) = g(x), x \in [0,1]"""
    return numpy.sin(math.pi*x)

#Initialize
if len(sys.argv) == 6:
    N       =   int(sys.argv[1])
    dt      = float(sys.argv[2])
    T       =   int(sys.argv[3])
    method  =   int(sys.argv[4])
    outfile =   str(sys.argv[5])
else:
    print "Usage:", sys.argv[0], "N dt T method"
    print "\t N:       Number of inner meshpoints"
    print "\t dt:      Length of timestep"
    print "\t T:       Number of timesteps"
    print "\t method:  Method selector:"
    print "\t          1 => Forward Euler"
    print "\t          2 => Backward Euler"
    print "\t          3 => Crank-Nicolson"
    print "\t outfile: Name of file to write to"
    sys.exit(0)

ofile = open(outfile,'w')

#dx = 1/float(N+1)
u = numpy.zeros((T,N+2),numpy.double)
(x,dx) = numpy.linspace (0,1,N+2, retstep=True)
alpha = dt/(dx**2)

#Initial codition
u[0,:] = g(x)
u[0,0] = u[0,N+1] = 0.0 #Implement boundaries rigidly

if   method == 1:
    forward_euler(alpha,u,N,T)
elif method == 2:
    backward_euler(alpha,u,N,T)
elif method == 3:
    crank_nicolson(alpha,u,N,T)
else:
    print "Please select method 1,2, or 3!"
    import sys
    sys.exit(0)

#Output
ostring = "%16.8E "*(N+2) + "\n";
for t in xrange(T):
    ofile.write(ostring % tuple(u[t,:]))
ofile.close()
