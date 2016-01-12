# coding=utf-8
#Brute force Monte Carlo example parallelized with pp
#Written by Kyrre Ness Sjøbæk

import math, numpy, sys
import pp

def func(x):
    """Integrand"""
    return 4/(1.0+x*x)

def montecarlo(N):
    mc_sum  = 0.0
    mc_sum2 = 0.0
    for i in xrange(N):
        fx=func(numpy.random.random())
        mc_sum  += fx
        mc_sum2 += fx*fx
    return (mc_sum, mc_sum2)

#Read in number of samples
if len(sys.argv) == 2:
    N = int(sys.argv[1])
else:
    print "Usage:",sys.argv[0],"N"
    sys.exit(0)

#Setup pp (only local CPU's)
job_server = pp.Server()

#Want number of samples to be divisible by number of cpus:
if N%job_server.get_ncpus() != 0:
    N += job_server.get_ncpus() -  N%job_server.get_ncpus()

print "Now running with",job_server.get_ncpus(),"nodes,",\
      N,"samples =>", N/job_server.get_ncpus(),"samples pr. node"

#Evaluate the integral using a crude MonteCarlo
#This is actually not the optimal way of using pp. pp can be run
#with more jobs than cpus to distribute the jobs more efficiently
#(see other parallel code in this chapter)
job = [];
for cpu in xrange(job_server.get_ncpus()):
    job.append(job_server.submit(montecarlo, args=(N/4,),\
                                 depfuncs=(func,),modules=("numpy",)))

#Get back stats
job_server.wait()
mc_sum  = 0.0
mc_sum2 = 0.0
for cpu in xrange(job_server.get_ncpus()):
    (mc_sum_local,mc_sum2_local) = job[cpu]()
    mc_sum  += mc_sum_local
    mc_sum2 += mc_sum2_local

#Fix the statistics
mc_sum  /= float(N)
mc_sum2 /= float(N)
variance = mc_sum2 - mc_sum**2

print "Integral=",mc_sum,", variance=",variance
job_server.print_stats()
