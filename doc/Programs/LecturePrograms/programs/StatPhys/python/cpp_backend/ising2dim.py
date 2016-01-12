import sys, numpy
import ising2dim_backend

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
temps = numpy.arange(temp_init,temp_end+temp_step/2,temp_step)

for temp in temps:
    answerlist = ising2dim_backend.monteCarlo(temp,size,trials)
    E_av       = answerlist[0]
    E_variance = answerlist[1]
    M_av       = answerlist[2]
    M_variance = answerlist[3]
    Mabs_av    = answerlist[4]
    #Use "tail -f <outfilename>" to se output real-time
    ofile.write("%15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n" % (temp, E_av, E_variance, M_av, M_variance, Mabs_av))
    
ofile.close()
