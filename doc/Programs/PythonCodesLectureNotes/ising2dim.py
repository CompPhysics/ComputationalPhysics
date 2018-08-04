from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import math, sys

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

def monteCarlo(temp, NSpins, MCcycles):
    """
    Calculate the energy and magnetization
    (\"straight\" and squared) for a given temperature

    Input:
    - temp:   Temperature to calculate for
    - NSpins:   dimension of square matrix
    - MCcycles: Monte-carlo MCcycles (how many times do we
                                  flip the matrix?)

    Output:
    - E_av:       Energy of matrix averaged over MCcycles, normalized to spins**2
    - E_variance: Variance of energy, same normalization * temp**2
    - M_av:       Magnetic field of matrix, averaged over MCcycles, normalized to spins**2
    - M_variance: Variance of magnetic field, same normalization * temp
    - Mabs:       Absolute value of magnetic field, averaged over MCcycles
    """

    #Setup spin matrix, initialize to ground state
    spin_matrix = np.zeros( (NSpins,NSpins), np.int8) + 1

    #Create and initialize variables
    E    = M     = 0
    E_av = E2_av = M_av = M2_av = Mabs_av = 0
    
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
                M += 2*spin_matrix[x,y]
                E += deltaE
            
        #Update expectation values
        E_av    += E
        E2_av   += E**2
        M_av    += M
        M2_av   += M**2
        Mabs_av += int(math.fabs(M))


    #Normalize average values
    E_av       /= float(MCcycles);
    E2_av      /= float(MCcycles);
    M_av       /= float(MCcycles);
    M2_av      /= float(MCcycles);
    Mabs_av    /= float(MCcycles);
    #Calculate variance and normalize to per-point and temp
    E_variance  = (E2_av-E_av*E_av)/float(NSpins*NSpins*temp*temp);
    M_variance  = (M2_av-M_av*M_av)/float(NSpins*NSpins*temp);
    #Normalize returned averages to per-point
    E_av       /= float(NSpins*NSpins);
    M_av       /= float(NSpins*NSpins);
    Mabs_av    /= float(NSpins*NSpins);
    
    return (E_av, E_variance, M_av, M_variance, Mabs_av)
    
    
# Main program
# temperature steps, initial temperature, final temperature
NumberTsteps = 20
InitialT = 1.5
FinalT = 2.5
Tsteps = (FinalT-InitialT)/NumberTsteps
Temp = np.zeros(NumberTsteps)
for T in range(NumberTsteps): 
    Temp[T] = InitialT+T*Tsteps
# Declare arrays that hold averages
Energy       = np.zeros(NumberTsteps);   Magnetization  = np.zeros(NumberTsteps)
SpecificHeat = np.zeros(NumberTsteps);   Susceptibility = np.zeros(NumberTsteps)
MagnetizationAbs = np.zeros(NumberTsteps)
# Define number of spins
NSpins = 20
# Define number of Monte Carlo cycles
MCcycles = 100000
# Perform the simulations over a range of temperatures
for T in range(NumberTsteps): 
    (Energy[T], SpecificHeat[T], Magnetization[T], Susceptibility[T], MagnetizationAbs[T]) = monteCarlo(Temp[T],NSpins,MCcycles)
# And finally plot
f = plt.figure(figsize=(18, 10)); # plot the calculated values    

sp =  f.add_subplot(2, 2, 1 );
plt.plot(Temp, Energy, 'o', color="green");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Energy ", fontsize=20);

sp =  f.add_subplot(2, 2, 2 );
plt.plot(Temp, abs(Magnetization), 'o', color="red");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Magnetization ", fontsize=20);

sp =  f.add_subplot(2, 2, 3 );
plt.plot(Temp, SpecificHeat, 'o', color="blue");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Specific Heat ", fontsize=20);

sp =  f.add_subplot(2, 2, 4 );
plt.plot(Temp, Susceptibility, 'o', color="black");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Susceptibility", fontsize=20);

plt.show()

