#Variational Monte Carlo program for the Helium atom which utilizes the code in MC.cpp - parallelized with pypar
#This implementation is based on a master-slave strategy for distributing work, instead of splitting up
#each integral on all the nodes.
#Written by Magnar K. Bugge

from math import sqrt
from numpy import linspace
import pypar
import MC

nprocs = pypar.size() #Number of processes      
myid = pypar.rank() #Id of this process

class result:
    alpha = 0.0
    E = 0.0
    sigma = 0.0
    error = 0.0
    acceptance = 0.0
    id = 0 #id of the process who did this job

MCcycles = 100000000 #Number of MC cycles
MCcycles2 = 10000 #Number of MC cycles for determination of optimal delta
delta_min = .01 #Minimum length of Metropolis step
delta_max = 2.0 #Maximum length of Metropolis step
tolerance = .01
idum = MC.seed() * (myid+1) #Seed for random number generator (different for each process)

#Function which should be close to zero for optimal delta
def difference(delta):
    x = MC.runMC(MCcycles2,delta,idum,alpha)
    return x.accepted*1.0/MCcycles2 - .5 #We want 50% accepted moves

#Function for sorting the results from small to large alpha
def sort(results):
    sorted_results = []
    
    for i in xrange(len(values)):
        sorted_results.append(result())
        sorted_results[i].alpha = values[i]
        j = 0
        while results[j].alpha != sorted_results[i].alpha:
            j += 1
        sorted_results[i].E = results[j].E
        sorted_results[i].sigma = results[j].sigma
        sorted_results[i].error = results[j].error
        sorted_results[i].acceptance = results[j].acceptance
        #add a line here if you also want to keep information about who (id) did what
    
    return sorted_results

#Array of alpha values
values = linspace(1.4,2.5,23) #(alpha values)

#Code for the master
if myid == 0:

    results = [] #List of results
    ndistributed = 0 #Number of distributed jobs
    nreceived = 0 #Number of received jobs

    #The master distributes jobs to the slaves. First, he sends one job to each slave.
    for i in xrange(1,nprocs):
        pypar.send(values[ndistributed],destination=i)
        ndistributed += 1

    #Then he goes in a loop, waiting for results and sending out new jobs
    while nreceived < len(values):

        #Receive a result
        results.append(pypar.receive(source=pypar.any_source))
        nreceived += 1

        #If there is more jobs to do, send the next one, if not, send a termination signal (alpha = 0)
        if ndistributed < len(values):
            pypar.send(values[ndistributed],destination=results[-1].id)
            ndistributed += 1
        else:
            pypar.send(0.0,destination=results[-1].id)

    #As the results are generally received out of order, they need to be sorted
    results = sort(results)

    #Finally, master writes the results to file
    outfile = open('data','w')
    for i in xrange(len(results)):
        outfile.write('%f %f %f %f %f\n' %(results[i].alpha,results[i].E,results[i].sigma,results[i].error,results[i].acceptance))
    outfile.close()

#Code for the slaves
if myid != 0:

    #Run in a loop receiving and returning jobs as long as the termination signal alpha=0 is not given
    while True:
        
        #Receive a job (or possibly the termination signal)
        alpha = pypar.receive(source=0)

        #If it isn't the termination signal, do the job and return the result
        if alpha != 0.0:

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
            x = MC.runMC(MCcycles,delta,idum,alpha)
            idum = x.idum

            #Calculate statistics and store them in a result object
            thisresult = result()
            thisresult.alpha = alpha
            thisresult.id = myid
            
            thisresult.E = x.sum / x.N
            E2 = x.squaresum / x.N
            thisresult.sigma = sqrt(E2 - thisresult.E**2)
            thisresult.acceptance = x.accepted*1.0/MCcycles
            thisresult.error = thisresult.sigma / sqrt(x.N)

            #Send to master
            pypar.send(thisresult,destination=0)
            
        else:
            break


pypar.finalize()
