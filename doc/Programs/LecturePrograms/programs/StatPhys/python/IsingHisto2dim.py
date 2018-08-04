from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import math, sys

def periodic (i, limit, add):
    return (i+limit+add) % limit

def monteCarlo(Energies, temp, NSpins, MCcycles):

    #Setup spin matrix, initialize to ground state
    spin_matrix = np.zeros( (NSpins,NSpins), np.int8) + 1

    E = M = 0.0    
    #Setup array for possible energy changes
    w = np.zeros(17,np.float64)
    for de in range(-8,9,4): #include +8
        w[de+8] = math.exp(-de/temp)
    
    #Calculate initial magnetization:
    M = spin_matrix.sum()
    #Calculate initial energy
    for j in range(NSpins): 
        for i in range(NSpins):
            E -= spin_matrix.item(i,j)*\
                 (spin_matrix.item(periodic(i,NSpins,-1),j) + spin_matrix.item(i,periodic(j,NSpins,1)))

    #Start metropolis MonteCarlo computation 
    for i in range(MCcycles):
        #Metropolis
        #Loop over all spins, pick a random spin each time
        for s in range(NSpins**2):
            x = int(np.random.random()*NSpins)
            y = int(np.random.random()*NSpins)
            deltaE = 2*spin_matrix.item(x,y)*\
                     (spin_matrix.item(periodic(x,NSpins,-1), y) +\
                      spin_matrix.item(periodic(x,NSpins,1),  y) +\
                      spin_matrix.item(x, periodic(y,NSpins,-1)) +\
                      spin_matrix.item(x, periodic(y,NSpins,1)))
            if np.random.random() <= w[deltaE+8]:
                #Accept!
                spin_matrix[x,y] *= -1
                E += deltaE
            
        #Update expectation values
        Energies[i]    += E

# Main program

# Define number of spins
NSpins = 20
# Define number of Monte Carlo cycles
MCcycles = 10000
# temperature steps, initial temperature, final temperature
Temp = 2.5
# Declare arrays that hold averages
Energies       = np.zeros(MCcycles)
# Obtain the energies to construct the diagram
monteCarlo(Energies,Temp,NSpins,MCcycles)

n, bins, patches = plt.hist(Energies, 100, facecolor='green')

plt.xlabel('$E$')
plt.ylabel('Energy distribution P(E)')
plt.title(r'Energy distribution at  $k_BT=2.5$')
plt.axis([-800, -300, 0, 500])
plt.grid(True)
plt.show()
