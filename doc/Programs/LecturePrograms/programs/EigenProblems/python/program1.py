#Program wich solves the one-particle (time independent) Schrodinger
#equation for a potential specified in function potential(). This
#example is for the hydrogen atom. This is a translation of the C++
#code for chapter 12 in "Computational Physics", M. H. Jensen
#Translated to Python by Magnar K. Bugge

import numpy
import math
import sys
import computationalLib

#Read output file name from command line argument, abort if none is given
if len(sys.argv) == 2:
    outfilename = sys.argv[1]
else:
    print '\nError: name of output file must be given as command line argument.\n'
    sys.exit(0)
    
#Function for initialisation of parameters
def initialize():
    r_min = eval(raw_input('Min value of r: '))
    r_max = eval(raw_input('Max value of r: '))
    orb_l = eval(raw_input('Orbital angular momentum (l): '))
    max_step = eval(raw_input('Number of steps: '))

    return r_min,r_max,orb_l,max_step

#Function for outputting results to file
def output():
    outfile = open(outfilename,'w')
    outfile.write('Results:\nr_min = %f\nr_max = %f\nOrbital angular momentum (l) = \
%d\nNumber of steps = %d\nFive lowest eigenvalues:\n' %(r_min,r_max,orb_l,max_step))
    for i in xrange(5):
        outfile.write('%f\n' %(d[i]))

#Potential function (returns value of potential for given r)
def potential(r):

    return -2.0/r


#Here starts the main program

#Read in data from user
r_min,r_max,orb_l,max_step = initialize()

#Initialize constants
step = (r_max - r_min) / max_step
const_2 = -1.0 / step**2
const_1 = -2.0 * const_2
orb_factor = orb_l * (orb_l + 1)

#Calculate array of potential values
V = numpy.zeros(max_step+1,numpy.double)
for i in xrange(max_step+1):
    r = r_min + i*step
    V[i] = potential(r) + orb_factor/r**2 #(include centrifugal term)

#Calculate elements of the matrix to be diagonalized
d = numpy.zeros(max_step,numpy.double) #diagonal elements
e = numpy.zeros(max_step,numpy.double) #off-diagonal elements
z = numpy.zeros((max_step,max_step),numpy.double) #matrix for eigenvectors (only used as a dummy-argument to
                                                  #tqli() in the current version of this program)
for i in xrange(max_step):
    d[i] = const_1 + V[i+1]
    e[i] = const_2

    z[i,i] = 1 #(identity matrix)

#Diagonalize and obtain eigenvalues. The eigenvalues are stored in d.
computationalLib.pylib(cpp=False).tqli(d,e,z)

#Sort the eigenvalues (smallest to largest)
d = numpy.msort(d)

#Output to file
output()
