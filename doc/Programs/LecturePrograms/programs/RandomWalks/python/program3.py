# coding=utf-8
#
# 1D-randomwalk: A walker makes several steps,
# with a given number of walks pr. trial.
# Similar to program1, but this version also cumputes
# the probability of a gridpoint being "touched"
#
#Translated to Python by Kyrre Ness Sjøbæk

import numpy, sys

def mc_trial(number_walks,move_probability,walk_cum,walk2_cum, probability):
    """
    Do a MonteCarlo trial, that is,
    random-walk one particle.

    Input:
    - number_walks:     Number of steps to walk the particle
    - move_probability: Probability that the particle
                        will step right when doing a step
    - walk_cum:         Numpy-array of length number_walks + 1,
                        containing the sum of the position
                        of the particles as a function of time
                        (usefull to calculate mean pos. as a function
                        of time)
    - walk2_cum:        Same as walk_cum, but with the sum of the
                        positions squared
    - probability:      Number of times each gridpoint is hit

    Output: As walk_cum, walk2_cum, and probability are (pointers to)
    numpy arrays, they are altered also in the calling function.
    """
    #Initial pos. As walk_cum[0]=walk2_cum[0] = 0.0
    #by initialization, it is uneccessary to add this step to
    #the arrays...
    pos = 0;
    
    for walk in xrange(number_walks+1):
        if numpy.random.random() <= move_probability:
            pos += 1
        else:
            pos -= 1
        walk_cum[walk]   += pos
        walk2_cum[walk]  += pos**2
        #Zero-position of the array is the leftmost
        #end of the grid
        probability[pos+number_walks] += 1

def mc_sample(trials, number_walks, move_probability):
    """
    Run as many trials as asked for by input variable trials.
    Wrapper to mc_trial, split out for easier paralellization

    Output:
    NumPy arrays walk_cum and walk2_cum, length number_walks + 1
    Numpy array probability, length 2*number_walks +1, which is
    equalient to the size of the reachable grid
    """

    walk_cum  = numpy.zeros(number_walks+1)
    walk2_cum = numpy.zeros(number_walks+1)
    
    probability = numpy.zeros(2*number_walks+1)
    
    for trial in xrange(trials):
        mc_trial(number_walks,move_probability,walk_cum,walk2_cum,probability)

    return (walk_cum,walk2_cum,probability)

#Main program

#Get data from command-line arguments
if len(sys.argv) == 6:
    outfilename          =       sys.argv[1]
    probability_filename =       sys.argv[2]
    trials               =   int(sys.argv[3])
    number_walks         =   int(sys.argv[4])
    move_probability     = float(sys.argv[5])
else:
    print "Usage: python", sys.argv[0],\
          "outfilename probability_filename trials number_walks move_probability"
    sys.exit(0)

#Open files at once, or fail loudly
ofile      = open(outfilename, 'w')
ofile_prob = open(probability_filename,'w')

#Do the MC
(walk_cum,walk2_cum,probability) = mc_sample(trials,number_walks, move_probability);

#Output
for i in xrange(len(walk_cum)):
    #Normalize to number of trials (= number of walkers)
    xaverage  = walk_cum[i]/float(trials)
    x2average = walk2_cum[i]/float(trials)
    variance = x2average - xaverage**2
    ofile.write("%6d %15.8E %15.8E\n" % (i,xaverage,variance))
ofile.close()

#Mean number of times each gridpoint is hit
norm = float(numpy.sum(probability))
for i in xrange(len(probability)):
    ofile_prob.write("%6d %15.8E\n" % (i+number_walks,probability[i]/norm))
ofile_prob.close()
