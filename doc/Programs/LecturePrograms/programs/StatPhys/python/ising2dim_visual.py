# coding=utf-8
#2-dimensional ising model with visualization
#Written by Kyrre Ness Sjøbæk

import numpy, sys, math
import pygame

#Needed for visualize when using SDL
screen = None;
font   = None;
BLOCKSIZE = 10

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

def visualize(spin_matrix, temp, E, M, method):
    """
    Visualize the spin matrix

    Methods:
    method = -1:No visualization (testing)
    method = 0: Just print it to the terminal
    method = 1: Pretty-print to terminal
    method = 2: SDL/pygame single-pixel
    method = 3: SDL/pygame rectangle
    """

    #Simple terminal dump
    if method == 0:
        print "temp:", temp, "E:", E, "M:", M
        print spin_matrix
    #Pretty-print to terminal
    elif method == 1:
        out = ""
        size = len(spin_matrix)
        for y in xrange(size):
            for x in xrange(size):
                if spin_matrix.item(x,y) == 1:
                    out += "X"
                else:
                    out += " "
            out += "\n"
        print "temp:", temp, "E:", E, "M:", M
        print out + "\n"
    #SDL single-pixel (useful for large arrays)
    elif method == 2:
        size = len(spin_matrix)
        screen.lock()
        for y in xrange(size):
            for x in xrange(size):
                if spin_matrix.item(x,y) == 1:
                    screen.set_at((x,y),(255,255,255))
                else:
                    screen.set_at((x,y),(0,0,0))
        screen.unlock()
        pygame.display.flip()
    #SDL block (usefull for smaller arrays)
    elif method == 3:
        size = len(spin_matrix)
        screen.lock()
        for y in xrange(size):
            for x in xrange(size):
                if spin_matrix.item(x,y) == 1:
                    rect = pygame.Rect(x*BLOCKSIZE,y*BLOCKSIZE,BLOCKSIZE,BLOCKSIZE)
                    pygame.draw.rect(screen,(255,255,255),rect)
                else:
                    rect = pygame.Rect(x*BLOCKSIZE,y*BLOCKSIZE,BLOCKSIZE,BLOCKSIZE)
                    pygame.draw.rect(screen,(0,0,0),rect)
        screen.unlock()
        pygame.display.flip()
    #SDL block w/ data-display
    elif method == 4:
        size = len(spin_matrix)
        screen.lock()
        for y in xrange(size):
            for x in xrange(size):
                if spin_matrix.item(x,y) == 1:
                    rect = pygame.Rect(x*BLOCKSIZE,y*BLOCKSIZE,BLOCKSIZE,BLOCKSIZE)
                    pygame.draw.rect(screen,(255,255,255),rect)
                else:
                    rect = pygame.Rect(x*BLOCKSIZE,y*BLOCKSIZE,BLOCKSIZE,BLOCKSIZE)
                    pygame.draw.rect(screen,(0,0,0),rect)
        s = font.render("<E> = %5.3E; <M> = %5.3E" % E,M,False,(255,0,0))
        screen.blit(s,(0,0))
        
        screen.unlock()
        pygame.display.flip()    

    

def monteCarlo(temp, size, trials, visual_method):
    """
    Calculate the energy and magnetization
    (\"straight\" and squared) for a given temperature

    Input:
    - temp:   Temperature to calculate for
    - size:   dimension of square matrix
    - trials: Monte-carlo trials (how many times do we
                                  flip the matrix?)
    - visual_method: What method should we use to visualize?

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

        visualize(spin_matrix, temp,E/float(size**2),M/float(size**2), method);

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
if len(sys.argv) == 5:
    size        =   int(sys.argv[1])
    trials      =   int(sys.argv[2])
    temp        = float(sys.argv[3])
    method      =   int(sys.argv[4])
else:
    print "Usage: python",sys.argv[0],\
          "lattice_size trials temp method"
    sys.exit(0)

if method > 4:
    print "method < 3!"
    sys.exit(0)

#Initialize pygame
if method == 2 or method == 3 or method == 4:
    pygame.init()
    if method == 2:
        screen = pygame.display.set_mode((size,size))
    elif method == 3:
        screen = pygame.display.set_mode((size*10,size*10))
    elif method == 4:
        screen = pygame.display.set_mode((size*10,size*10))
        font   = pygame.font.Font(None,12)
        

(E_av, E_variance, M_av, M_variance, Mabs_av) = monteCarlo(temp,size,trials, method)
print "%15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n" % (temp, E_av, E_variance, M_av, M_variance, Mabs_av)

pygame.quit();
