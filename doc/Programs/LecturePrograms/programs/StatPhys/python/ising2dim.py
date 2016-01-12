# coding=utf-8
#2-dimensional Ising model
#Translated to Python by Kyrre Ness Sjøbæk

import numpy, sys, math


def periodic (i, limit, add):
    """
    Choose correct matrix index with periodic
    boundary conditions

    Input:
    - i:     Base index
    - limit: Highest \"legal\" index
    - add:   Number to add or subtract from i
    """
    return (i+limit+add) % limit

def monteCarlo(temp, size, trials):
    """
    Calculate the energy and magnetization
    (\"straight\" and squared) for a given temperature

    Input:
    - temp:   Temperature to calculate for
    - size:   dimension of square matrix
    - trials: Monte-carlo trials (how many times do we
                                  flip the matrix?)

    Output:
    - E_av:       Energy of matrix averaged over trials, normalized to spins**2
    - E_variance: Variance of energy, same normalization * temp**2
    - M_av:       Magnetic field of matrix, averaged over trials, normalized to spins**2
    - M_variance: Variance of magnetic field, same normalization * temp
    - Mabs:       Absolute value of magnetic field, averaged over trials
    """

    #Setup spin matrix, initialize to ground state
    spin_matrix = numpy.zeros( (size,size), numpy.int8) + 1

    #Create and initialize variables
    E    = M     = 0
    E_av = E2_av = M_av = M2_av = Mabs_av = 0
    
    #Setup array for possible energy changes
    w = numpy.zeros(17,numpy.float64)
    for de in xrange(-8,9,4): #include +8
        w[de+8] = math.exp(-de/temp)
    
    #Calculate initial magnetization:
    M = spin_matrix.sum()
    #Calculate initial energy
    for j in xrange(size): 
        for i in xrange(size):
            E -= spin_matrix.item(i,j)*\
                 (spin_matrix.item(periodic(i,size,-1),j) + spin_matrix.item(i,periodic(j,size,1)))

    #Start metropolis MonteCarlo computation 
    for i in xrange(trials):
        #Metropolis
        #Loop over all spins, pick a random spin each time
        for s in xrange(size**2):
            x = int(numpy.random.random()*size)
            y = int(numpy.random.random()*size)
            deltaE = 2*spin_matrix.item(x,y)*\
                     (spin_matrix.item(periodic(x,size,-1), y) +\
                      spin_matrix.item(periodic(x,size,1),  y) +\
                      spin_matrix.item(x, periodic(y,size,-1)) +\
                      spin_matrix.item(x, periodic(y,size,1)))
            if numpy.random.random() <= w[deltaE+8]:
                #Accept!
                spin_matrix[x,y] *= -1
                M += 2*spin_matrix[x,y]
                E += deltaE
            
        #Update expectation values
        E_av    += E
        E2_av   += E**2
        M_av    += M
        M2_av   += M**2
        Mabs_av += int(math.fabs(M))

        #DEBUG
#         Etest = 0
#         for j in xrange(size): 
#             for i in xrange(size):
#                 Etest -= spin_matrix.item(i,j)*\
#                          (spin_matrix.item(periodic(i,size,-1),j) + spin_matrix.item(i,periodic(j,size,-1)))
#         print Etest
#         print E
#         print spin_matrix

    #Output results, normalized to number of spins
#     print "Temp:      ", temp
#     print "Trials:    ", trials
#     print "E_av:      ", E_av/float(trials*size**2)
#     print "E_variance:", (E2_av/float(trials) - E_av**2/float(trials**2))/float(size**2 * temp**2)
#     print "M_av:      ", M_av/float(trials*size**2)
#     print "M_variance:", (M2_av/float(trials) - M_av**2/float(trials**2))/float(size**2 * temp)
#     print "Mabs_av:   ", Mabs_av/float(trials*size**2)
#     print

    #Normalize average values
    E_av       /= float(trials);
    E2_av      /= float(trials);
    M_av       /= float(trials);
    M2_av      /= float(trials);
    Mabs_av    /= float(trials);
    #Calculate variance and normalize to per-point and temp
    E_variance  = (E2_av-E_av*E_av)/float(size*size*temp*temp);
    M_variance  = (M2_av-M_av*M_av)/float(size*size*temp);
    #Normalize returned averages to per-point
    E_av       /= float(size*size);
    M_av       /= float(size*size);
    Mabs_av    /= float(size*size);
    
    return (E_av, E_variance, M_av, M_variance, Mabs_av)
    
    
# Main program

#Get input
if len(sys.argv) == 7:
    outfilename =       sys.argv[1]
    size        =   int(sys.argv[2])
    trials      =   int(sys.argv[3])
    temp_init   = float(sys.argv[4])
    temp_end    = float(sys.argv[5])
    temp_step   = float(sys.argv[6])
else:
    print "Usage: python",sys.argv[0],\
          "outfilename lattice_size trials temp_init temp_end temp_step"

    sys.exit(0)

ofile = open(outfilename,'w')

#Loop over temperatures (highly advantagous to paralellize!)

#arange has round-off problems, and sec comes from scipy,
#which is hard to compile. Work-around-warning!
temps = numpy.arange(temp_init,temp_end+temp_step/2,temp_step,float)

for temp in temps:
    (E_av, E_variance, M_av, M_variance, Mabs_av) = monteCarlo(temp,size,trials)
    #Use "tail -f <outfilename>" to se output real-time
    ofile.write("%15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n" % (temp, E_av, E_variance, M_av, M_variance, Mabs_av))
    
ofile.close()
