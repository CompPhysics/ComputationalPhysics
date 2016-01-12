# coding=utf-8
#Radioactive decay of nuclei - parallel version
#Written by Kyrre Ness Sjøbæk

import sys, numpy
import pp

def montecarlo(num_particles_init, max_time, decay_prob):
    """Loop over max_time timesteps to generate one MC sample,
    return an array containing number of left particles as a
    function of timesteps"""
    ncum     = numpy.zeros(max_time+1,numpy.int)
    ncum[0]  = num_particles_init
    n_unstab = num_particles_init
    
    for t in xrange(1,max_time):
        for particle in xrange(n_unstab):
            if numpy.random.random() <= decay_prob:
                n_unstab -= 1
        ncum[t] = n_unstab
        if n_unstab == 0:
            break
    return ncum
    
#Read in inputdata
if len(sys.argv) == 6:
    outfilename        = sys.argv[1]
    num_particles_init = int(sys.argv[2])
    max_time           = int(sys.argv[3])
    num_cycles         = int(sys.argv[4])
    decay_prob         = float(sys.argv[5])
else:
    print "Usage:",sys.argv[0],"outfilename num_particles_init",\
          "max_time num_cycles decay_prob"
    sys.exit(0)

#Open the file before we do any calculations, in case this fails...
ofile = open(outfilename,'w')

#Setup pp
job_server=pp.Server()

#Dispatch jobs, one job = one MC-cycle
jobs = [job_server.submit(montecarlo, args=(num_particles_init,max_time,decay_prob),\
                          modules=("numpy",)) for cycle in xrange(num_cycles)]
#Get back results
ncum = numpy.zeros(max_time+1,numpy.int)
for job in jobs:
    ncum += job()

#Write results to file, normalized to number_particles
for i in xrange(len(ncum)):
    ofile.write("%E\n" % (ncum[i]/float(num_cycles)))

ofile.close()
