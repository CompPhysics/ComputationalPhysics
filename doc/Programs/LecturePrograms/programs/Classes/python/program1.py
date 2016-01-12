# -*- coding: utf-8 -*-

#This program calculates the 2nd derivative of exp(x)
#at a given x, for a range of stepsizes,
#using the three-point-formula
#Translated to Python by Kyrre Ness Sjøbæk

#Load libraries
from math import * #Dump everything in module "math" into the local namespace;
                   #To se what this means, import and use dir() in an interactive session

#Functions called by main program:
#Gets the initial stepsize h, 
def initialize():

    #We can import modules wherever we want!
    import sys; #System functions, such as reading from standard input
                #(usually keyboard)
    
    print "Initial stepsize:"
    initial_step    = float(sys.stdin.readline())

    print "Evaluate at point x:"
    x               = float(sys.stdin.readline())

    print "Number of steps (stepsize will be halved each iteration):"
    number_of_steps = int(sys.stdin.readline())

    #Return a tuple containing the interesting data
    return (initial_step, x, number_of_steps)

#Computes the second derivative of exp(x)
def second_derivative(number_of_steps, x, initial_step):
    h                   = initial_step
    h_step              = [] #Create an empty python list
    computed_derivative = []

    for i in xrange(number_of_steps):
        h_step.append(h)
        computed_derivative.append( (exp(x+h)-2*exp(x)+exp(x-h)) / (h*h) )
        h = h/2.0

    return (h_step,computed_derivative)

#Writes log10 of the stepsize and relative error to a file named "out.dat"
def output (h_step, computed_derivative, x):
    
    output_file = open("out.dat", 'w') #Open file for writing

    for i in xrange(len(h_step)):
        #C-like output formatting may be used -
        #Usefull when controll over num. format needed
        output_file.write('%g %12.5E\n' % ( log10(h_step[i]),\
                                            log10( (fabs(computed_derivative[i]-exp(x))) /exp(x)) ) )

    #Close the file, or we might not get all the data written...
    output_file.close()


#Main code, ran when program is started
#This has to be below the functions called in it!
(initial_step, x, number_of_steps) = initialize()
(h_step, computed_derivative)      = second_derivative(number_of_steps, x, initial_step)
output (h_step, computed_derivative, x)
