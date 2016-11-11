#Variational Monte Carlo program for the Helium atom - parallelized with pypar
#Written by Magnar K. Bugge

from math import sqrt,exp
from numpy import zeros,double,linspace,array
from random import random
import pypar

MCcycles = 10000000 #Number of MC cycles
MCcycles2 = 10000 #Number of MC cycles for determination of optimal delta
delta_min = .01 #Minimum length of Metropolis step
delta_max = 2.0 #Maximum length of Metropolis step
tolerance = .01

nprocs = pypar.size() #Number of processes      
myid = pypar.rank() #Id of this process

#Function which checks for singularity
def hasSingularity(R):
    r1 = sqrt(R[0,0]**2 + R[0,1]**2 + R[0,2]**2)
    r2 = sqrt(R[1,0]**2 + R[1,1]**2 + R[1,2]**2)
    r12 = sqrt((R[0,0]-R[1,0])**2 + (R[0,1]-R[1,1])**2 + (R[0,2]-R[1,2])**2)

    if (r1 < 1e-10 or r2 < 1e-10 or r12 < 1e-10):
        return True
    else:
        return False

#Local energy
def E_local(R):
    r1 = sqrt(R[0,0]**2 + R[0,1]**2 + R[0,2]**2)
    r2 = sqrt(R[1,0]**2 + R[1,1]**2 + R[1,2]**2)
    r12 = sqrt((R[0,0]-R[1,0])**2 + (R[0,1]-R[1,1])**2 + (R[0,2]-R[1,2])**2)

    return (alpha - 2)*(1/r1 + 1/r2) + 1/r12 - alpha**2

#Trial wave function
def Psi_trial(R):
    r1 = sqrt(R[0,0]**2 + R[0,1]**2 + R[0,2]**2)
    r2 = sqrt(R[1,0]**2 + R[1,1]**2 + R[1,2]**2)

    return exp(-alpha*(r1+r2))

#Monte Carlo function
def runMC(MCcycles,delta):
    sum = 0.0
    squaresum = 0.0
    N = 0
    accepted = 0

    #Initialize positions
    R = zeros((2,3),double)
    for i in xrange(2):
        for j in xrange(3):
            R[i,j] = .5 * (random()*2 - 1)

    R_trial = R.copy()

    #Loop over cycles
    for k in xrange(MCcycles):
        #Trial position
        for i in xrange(2):
            for j in xrange(3):
                R_trial[i,j] = R[i,j] + delta * (random()*2 - 1)

        P_trial = Psi_trial(R_trial)
        P_trial = P_trial**2
        P = Psi_trial(R)
        P = P**2

        #Metropolis test
        if P_trial > P:
            accepted += 1
            R = R_trial.copy()
        elif random() < P_trial / P:
            accepted += 1
            R = R_trial.copy()
            
        #Calculate contribution if we do not have singularity and are finished thermalising
        if (not hasSingularity(R)) and k > MCcycles/20:
            EL = E_local(R)
                
            sum += EL
            squaresum += EL**2
            N += 1
            
    return sum,squaresum,N,accepted
            
#Function which should be close to zero for optimal delta
def difference(delta):
    sum,squaresum,N,accepted = runMC(MCcycles2,delta)
    return accepted*1.0/MCcycles2 - .5 #We want 50% accepted moves

#Array of alpha values
values = linspace(1.4,2.0,13) #(alpha values)

if myid == 0:
    outfile = open('data','w')

#Loop over alpha values
for alpha in values:

    #Determination of optimal delta value (for each alpha), i.e.
    #finding the zero-point of the difference function by the bisection method
    minimum = delta_min
    maximum = delta_max
    while maximum - minimum > tolerance:
        if difference(minimum)*difference((minimum+maximum)/2) < 0:
            maximum = (minimum + maximum) / 2
        else:
            minimum = (minimum + maximum) / 2
    delta = (minimum + maximum) / 2
    
    #Run MC calculation
    sum,squaresum,N,accepted = runMC(MCcycles/nprocs,delta)

    #All nodes send their stats to the master
    if myid != 0:
        pypar.send((sum,squaresum,N,accepted),destination=0)

    #The master receives, calculates and writes to file
    if myid == 0:
        for i in xrange(1,nprocs):
            isum,isquaresum,iN,iaccepted = pypar.receive(source=i)
            sum += isum
            squaresum += isquaresum
            N += iN
            accepted += iaccepted

        #Calculate statistics
        E = sum / N
        E2 = squaresum / N
        sigma = sqrt(E2 - E**2)
        acceptance = accepted*1.0/MCcycles
        error = sigma / sqrt(N)

        #Print results to screen
        #print 'alpha = %f, <E> = %f, sigma = %f, acceptance = %f, error = %f' %(alpha,E,sigma,error,acceptance)

        #Print results to file
        outfile.write('%f %f %f %f %f\n' %(alpha,E,sigma,error,acceptance))

if myid == 0:
    outfile.close()

pypar.finalize()
