#Variational Monte Carlo program for the Helium atom which utilizes the code in MC.cpp - parallelized with pypar
#Written by Magnar K. Bugge

from math import sqrt
from numpy import linspace
import pypar
import MC

nprocs = pypar.size() #Number of processes      
myid = pypar.rank() #Id of this process

MCcycles = 10000000 #Number of MC cycles
MCcycles2 = 10000 #Number of MC cycles for determination of optimal delta
delta_min = .01 #Minimum length of Metropolis step
delta_max = 2.0 #Maximum length of Metropolis step
tolerance = .01
idum = MC.seed() * (myid+1) #Seed for random number generator (different for each process)

#Function which should be close to zero for optimal delta
def difference(delta):
    x = MC.runMC(MCcycles2,delta,idum,alpha)
    return x.accepted*1.0/MCcycles2 - .5 #We want 50% accepted moves

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
    
    #Run MC calculation (store results in x)
    x = MC.runMC(MCcycles/nprocs,delta,idum,alpha)
    idum = x.idum

    #All nodes send their stats to the master
    if myid != 0:
        pypar.send((x.sum,x.squaresum,x.N,x.accepted),destination=0)

    #The master receives, calculates and writes to file
    if myid == 0:
        for i in xrange(1,nprocs):
            isum,isquaresum,iN,iaccepted = pypar.receive(source=i)
            x.sum += isum
            x.squaresum += isquaresum
            x.N += iN
            x.accepted += iaccepted

        #Calculate statistics
        E = x.sum / x.N
        E2 = x.squaresum / x.N
        sigma = sqrt(E2 - E**2)
        acceptance = x.accepted*1.0/MCcycles
        error = sigma / sqrt(x.N)

        #Print results to screen
        #print 'alpha = %f, <E> = %f, sigma = %f, error = %f, acceptance = %f' %(alpha,E,sigma,error,acceptance)

        #Print results to file
        outfile.write('%f %f %f %f %f\n' %(alpha,E,sigma,error,acceptance))

if myid == 0:
    outfile.close()

pypar.finalize()
