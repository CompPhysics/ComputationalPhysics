# coding=utf-8
#
# 1D-randomwalk: A walker makes several steps,
# with a given number of walks pr. trial
#
#Translated to Python by Kyrre Ness Sjøbæk

import numpy, sys

def mc_trial(number_walks,move_probability,walk_cum,walk2_cum):
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

    Output: As walk_cum and walk2_cum are numpy arrays, they are altered.
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
        walk_cum[walk]  += pos
        walk2_cum[walk] += pos**2

def mc_sample(trials, number_walks, move_probability):
    """
    Run as many trials as asked for by input variable trials.
    Wrapper to mc_trial, split out for easier paralellization

    Output: NumPy arrays walk_cum and walk2_cum, length number_walks + 1
    """

    walk_cum  = numpy.zeros(number_walks+1)
    walk2_cum = numpy.zeros(number_walks+1)
    for trial in xrange(trials):
        mc_trial(number_walks,move_probability,walk_cum,walk2_cum)

    return (walk_cum,walk2_cum)

#Main program

#Get data from command-line arguments
if len(sys.argv) == 5:
    outfilename      =       sys.argv[1]
    trials           =   int(sys.argv[2])
    number_walks     =   int(sys.argv[3])
    move_probability = float(sys.argv[4])
else:
    print "Usage: python", sys.argv[0], "outfilename trials number_walks move_probability"
    sys.exit(0)

#Open file at once, or fail loudly
ofile = open(outfilename, 'w')

#Do the MC
(walk_cum,walk2_cum) = mc_sample(trials,number_walks, move_probability);

#Output
for i in xrange(len(walk_cum)):
    #Normalize to number of trials (= number of walkers)
    xaverage  = walk_cum[i]/float(trials)
    x2average = walk2_cum[i]/float(trials)
    variance = x2average - xaverage**2
    ofile.write("%6d %15.8E %15.8E\n" % (i,xaverage,variance))

ofile.close()
    
