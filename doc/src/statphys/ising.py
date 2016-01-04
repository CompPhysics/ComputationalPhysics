
import numpy, sys, math
from  matplotlib import pyplot as plt
import numpy as np


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
    """

    #Setup spin matrix, initialize to ground state
    spin_matrix = numpy.zeros( (size,size), numpy.int8) + 1

    #Create and initialize variables
    E = 0
    E_av = E2_av = 0
    
    #Setup array for possible energy changes
    w = numpy.zeros(17,numpy.float64)
    for de in xrange(-8,9,4): #include +8
        w[de+8] = math.exp(-de/temp)
    
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
                E += deltaE
            
        #Update expectation values
        E_av    += E
        E2_av   += E**2

    E_av       /= float(trials);
    E2_av      /= float(trials);
    #Calculate variance and normalize to per-point and temp
    E_variance  = (E2_av-E_av*E_av)/float(size*size*temp*temp);
    #Normalize returned averages to per-point
    E_av       /= float(size*size);

    return (E_av, E_variance)
    
    
# Main program

# values of the lattice, number of Monte Carlo cycles and temperature domain
size        =   10
trials      =   10000
temp_init   = 1.8
temp_end    = 2.6
temp_step   = 0.1


temps = numpy.arange(temp_init,temp_end+temp_step/2,temp_step,float)
Dim = np.size(temps)
energy = np.zeros(Dim)
heatcapacity = np.zeros(Dim) 
temperature = np.zeros(Dim)
for temp in temps:
    (E_av, E_variance) = monteCarlo(temp,size,trials)
    temperature[temp] = temp
    energy[temp] = E_av
    heatcapacity[temp] = E_variance
plt.figure(1)
plt.subplot(211)
plt.axis([1.8,2.6,-2.0, -1.0])
plt.xlabel(r'Temperature $J/(k_B)$')
plt.ylabel(r'Average energy per spin  $E/N$')
plt.plot(temperature, energy, 'b-')
plt.subplot(212)
plt.axis([1.8,2.6, 0.0, 2.0])
plt.plot(temperature, heatcapacity, 'r-')
plt.xlabel(r'Temperature $J/(k_B)$')
plt.ylabel(r'Heat capacity per spin  $C_V/N$')
plt.savefig('energycv.pdf')
plt.show()

