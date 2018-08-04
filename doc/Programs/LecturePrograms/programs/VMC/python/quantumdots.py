# Importing various packages
from math import exp, sqrt
from random import random, seed
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


#Trial wave function for quantum dots in two dims
def WaveFunction(r,alpha,beta):
    argument = 0.0
    wf = 0.0
    # standard harmonic oscillator
    for i in range(NumberParticles):
        r_single_particle = 0.0
        for j in range(Dimension):
            r_single_particle += r[i,j]**2
        argument += (r_single_particle)
    wf = exp(-0.5*argument*alpha)
    #Electron-electron contribution
    for i1 in range(NumberParticles-1):
        for i2 in range(i1+1,NumberParticles):
            r_12 = 0.0
            for j in range(Dimension):
                r_12 += (r[i1,j] - r[i2,j])**2
            wf *= exp(r_12/(1.0+beta*r_12))
    return wf

#Local energy  for quantum dots in two dims
def LocalEnergy(r,wf,alpha,beta):
    
 #Kinetic energy
    r_plus = r.copy()
    r_minus = r.copy()
    e_kinetic = 0.0
    for i in range(NumberParticles):
        for j in range(Dimension):
            r_plus[i,j] = r[i,j] + h
            r_minus[i,j] = r[i,j] - h
            wf_minus = WaveFunction(r_minus,alpha,beta)
            wf_plus = WaveFunction(r_plus,alpha,beta)
            e_kinetic -= wf_minus+wf_plus-2*wf;
            r_plus[i,j] = r[i,j]
            r_minus[i,j] = r[i,j]

    e_kinetic = .5*h2*e_kinetic/wf
    
    #Potential energy
    e_potential = 0.0
    
    #Harmonic oscillator contribution
    for i in range(NumberParticles):
        r_single_particle = 0.0
        for j in range(Dimension):
            r_single_particle += r[i,j]**2
        e_potential += 0.5*r_single_particle

    #Electron-electron contribution
    for i1 in range(NumberParticles-1):
        for i2 in range(i1+1,NumberParticles):
            r_12 = 0.0
            for j in range(Dimension):
                r_12 += (r[i1,j] - r[i2,j])**2
            e_potential += 1.0/sqrt(r_12)
    
    return e_potential + e_kinetic

# The Monte Carlo sampling with the Metropolis algo
def MonteCarloSampling():

    NumberMCcycles= 10000
    StepSize = 1.0
    # positions
    PositionOld = np.zeros((NumberParticles,Dimension), np.double)
    PositionNew = np.zeros((NumberParticles,Dimension), np.double)
    # seed for rng generator
    seed()
    # start variational parameter
    alpha = 0.9
    for ia in range(MaxVariations):
        alpha += .025
        AlphaValues[ia] = alpha
        beta = 0.3 
        for jb in range(MaxVariations):
            beta += .05
            BetaValues[jb] = beta
            energy = energy2 = 0.0
            DeltaE = 0.0
            #Initial position
            for i in range(NumberParticles):
                for j in range(Dimension):
                    PositionOld[i,j] = StepSize * (random() - .5)
            wfold = WaveFunction(PositionOld,alpha,beta)

            #Loop over MC MCcycles
            for MCcycle in range(NumberMCcycles):
                #Trial position
                for i in range(NumberParticles):
                    for j in range(Dimension):
                        PositionNew[i,j] = PositionOld[i,j] + StepSize * (random() - .5)
                wfnew = WaveFunction(PositionNew,alpha,beta)

                #Metropolis test to see whether we accept the move
                if random() < wfnew**2 / wfold**2:
                   PositionOld = PositionNew.copy()
                   wfold = wfnew
                   DeltaE = LocalEnergy(PositionOld,wfold,alpha,beta)
                energy += DeltaE
                energy2 += DeltaE**2

            #We calculate mean, variance and error ...
            energy /= NumberMCcycles
            energy2 /= NumberMCcycles
            variance = energy2 - energy**2
            error = sqrt(variance/NumberMCcycles)
            Energies[ia][jb] = energy    
    return Energies, AlphaValues, BetaValues
#         print(energy, alpha, beta)        


#Here starts the main program with variable declarations
h = 0.001
h2 = 1.0/(h*h)
NumberParticles = 2
Dimension = 2
MaxVariations = 10
Energies = np.zeros((MaxVariations,MaxVariations))
AlphaValues = np.zeros(MaxVariations)
BetaValues = np.zeros(MaxVariations)
(Energies, AlphaValues, BetaValues) = MonteCarloSampling()
fig = plt.figure()
ax = fig.gca(projection='3d')
# Plot the surface.
X, Y = np.meshgrid(AlphaValues, BetaValues)
surf = ax.plot_surface(X, Y, Energies,cmap=cm.coolwarm,linewidth=0, antialiased=False)
# Customize the z axis.
ax.set_zlim(3.0, 3.2)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()


