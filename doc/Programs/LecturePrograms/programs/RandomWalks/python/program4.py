# coding=utf-8
#
# 1D-randomwalk: A walker makes several steps,
# with a given number of walks pr. trial.
# Similar to program1, but this version also cumputes
# the probability of a gridpoint being "touched"
#
#Translated to Python by Kyrre Ness Sjøbæk

import numpy, sys, math

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

def mc_sample(length,trials, number_walks, move_probability):
    """
    Generate the probability distribution for finding
    a walker at a gridpoint, after a number of walks on a
    1d lattice with wrap-around boundary conditions
    (\"pacman universe\")

    Input:
    - length: Lattice-points away from x=0
    - trials: Number of MonteCarlo trials (number of walkers)
    - move_probability: Probability of moving right

    Output:
    Normalized probability of finding a walker on a
    specific grid position
    """

    #Grid position of every walker
    x = numpy.zeros(trials,numpy.int)

    #Loop over timesteps and walkers,
    #and find the walkers "ending positions"
    for t in xrange(number_walks):
        for i in xrange(trials):
            if numpy.random.random() <= move_probability:
                x[i] += 1
                #Wraparound?
                if x[i] > length:
                    x[i] = -length
            else:
                x[i] -= 1
                if x[i] < -length:
                    x[i] = +length
    
    #Calculate the probability of finding a walker
    #each grid-position
    probability = numpy.zeros(2*length+1)
    for i in xrange(len(probability)):
        pos = i-length
        #count number of occurences of this pos i x array
        count = 0
        for j in xrange(len(x)):
            if x[j] == pos:
                count += 1
        #Normalize and save
        probability[i] = count/float(trials)
    return probability

#Main program

#Get data from command-line arguments
if len(sys.argv) == 5:
    length               =   int(sys.argv[1])
    trials               =   int(sys.argv[2])
    number_walks         =   int(sys.argv[3])
    move_probability     = float(sys.argv[4])
else:
    print "Usage: python", sys.argv[0],\
          "lenght trials number_walks move_probability"
    sys.exit(0)

#Do the MC
probability = mc_sample(length,trials,number_walks,move_probability);

#Not reliable: ln(0)
#entropy = - numpy.sum(probability*numpy.log(probability))

entropy = 0.0
for i in xrange(len(probability)):
    if probability[i] > 0.0:
        entropy -= probability[i]*math.log(probability[i])

print "Timesteps             =",number_walks
print "Walkers (num. trials) =",trials
print "Entropy               =",entropy
print
if len(probability) <= 101:
    print "Probability distribution (Flat => high entropy):"
    print probability
else:
    print "Probability distribution to big to print"
