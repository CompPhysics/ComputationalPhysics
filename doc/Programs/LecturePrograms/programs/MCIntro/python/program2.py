# coding=utf-8
#Particles in box example
#Translated to Python by Kyrre Ness Sjøbæk

import sys, numpy.random

#Get name of output file and number of particles in system
if len(sys.argv) == 3:
    filename = sys.argv[1]
    N        = int(sys.argv[2])
else:
    print "Usage:",sys.argv[0],"filename N"
    sys.exit(0)

ofile = open(filename, 'w')

#All particles start in left-hand box
Nl = N

#Do the MC!
for t in xrange(10*N):
    if (int(N*numpy.random.random()) <= Nl):
        Nl -= 1
    else:
        Nl += 1
    ofile.write("%g %g\n" % (t,Nl))

ofile.close()
