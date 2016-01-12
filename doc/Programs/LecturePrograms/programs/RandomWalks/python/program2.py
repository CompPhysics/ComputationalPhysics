# coding=utf-8
#
# 2D-randomwalk: A walker makes several steps,
# with a given number of walks pr. trial
#
#Translated to Python by Kyrre Ness Sjøbæk

import numpy, sys
#import math

def mc_trial(number_walks,move_probability,\
             walk_cum_x,walk_cum_y,walk2_cum_x,walk2_cum_y):
    """
    Do a MonteCarlo trial, that is,
    random-walk one particle.

    Input:
    - number_walks:     Number of steps to walk the particle
    - move_probability: Probability that the particle
                        will step right when doing a step.
                        See code for how other probs are
                        computed as a function of this.
    - walk_cum_x:       Numpy-array of length number_walks + 1,
                        containing the sum of the position
                        of the particles as a function of time
                        (usefull to calculate mean pos. as a function
                        of time)
    - walk_cum_y:       --**-- but for y
    - walk2_cum_x:      Same as walk_cum, but with the sum of the
                        positions squared
    - walk2_cum_y:      --**-- but for y

    Output: As walk_cum and walk2_cum are numpy arrays, they are altered.
    """
    #Initial pos. As walk_cum[0]=walk2_cum[0] = 0.0
    #by initialization, it is uneccessary to add this step to
    #the arrays...
    pos_x = 0;
    pos_y = 0;
    
    for walk in xrange(number_walks+1):
        rand = numpy.random.random();
        if rand <= move_probability:
            pos_x += 1
        elif rand <= 2*move_probability:
            pos_x -= 1
        elif rand <= 3*move_probability:
            pos_y += 1
        else:
            pos_y -= 1
        walk_cum_x[walk]  += pos_x
        walk_cum_y[walk]  += pos_y
        walk2_cum_x[walk] += pos_x**2
        walk2_cum_y[walk] += pos_y**2

def mc_sample(trials, number_walks, move_probability):
    """
    Run as many trials as asked for by input variable trials.
    Wrapper to mc_trial, split out for easier paralellization

    Output: NumPy arrays walk_cum_x,walk_cum_y,
    walk2_cum_x, and walk2_cum_y, length number_walks + 1
    """

    walk_cum_x  = numpy.zeros(number_walks+1)
    walk_cum_y  = numpy.zeros(number_walks+1)
    walk2_cum_x = numpy.zeros(number_walks+1)
    walk2_cum_y = numpy.zeros(number_walks+1)
    for trial in xrange(trials):
        mc_trial(number_walks,move_probability,\
                 walk_cum_x,walk_cum_y,walk2_cum_x,walk2_cum_y)

    return (walk_cum_x,walk_cum_y,walk2_cum_x,walk2_cum_y)

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
(walk_cum_x,walk_cum_y,walk2_cum_x,walk2_cum_y) = mc_sample(trials,number_walks, move_probability);

#Output
for i in xrange(len(walk_cum_x)):
    #Normalize to number of trials (= number of walkers)
    xaverage  = walk_cum_x[i]/float(trials)
    x2average = walk2_cum_x[i]/float(trials)
    xvariance = x2average - xaverage**2

    yaverage  = walk_cum_y[i]/float(trials)
    y2average = walk2_cum_y[i]/float(trials)
    yvariance = y2average - yaverage**2

    #total_average  = math.sqrt(xaverage**2 + yaverage**2)
    total_average  = xaverage + yaverage
    total_variance = xvariance + yvariance
    
    ofile.write("%6d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n" % (i,xaverage,xvariance,yaverage,yvariance,total_average,total_variance))

ofile.close()
    
