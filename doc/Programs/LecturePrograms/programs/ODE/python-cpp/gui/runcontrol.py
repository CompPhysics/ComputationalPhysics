#Do stuff so we can import ../swig/pendelum.py
import os,sys
sys.path.append(os.path.join(os.path.abspath(".."), "swig"))
from pendelum import pendelum

import shutil

debug = True

class runset:
    """
    This class represents a set of runs,
    which is a fixed set of (simulation-) parameters
    plus info about which methods have been run on these sets.

    Note that N is set per-run.
    """

    #Our list of runs. Key is method-runname
    runs = {}

    def __init__(self, basedir, runsetname,
                 m,l,omega,A,viscosity,phi_0,v_0,t_end):
        """
        Constructor, this gets the basic parameters
        for the simulations, as well as the name of the run,
        and where to store the dataset(s).

        The dataset is stored in a folder of the same name as
        runsetname, created under basedir.

        Input: (the non-obvious ones):
         - basedir:    Where should our runset folder be created
         - runsetname: What is the name of our runset
         - t_end:      Final time AS A MULTIPLUM OF PI
        """
        self.m         = float(m)
        self.l         = float(l)
        self.omega     = float(omega)
        self.A         = float(A)
        self.viscosity = float(viscosity)
        self.phi_0     = float(phi_0)
        self.v_0       = float(v_0)
        self.t_end     = float(t_end)

        #Setup paths etc., sanity checks
        if not os.path.isdir(basedir):
            print "Oops! Basedir doesn't exist!"
            sys.exit(0) #TODO: Should prob. use some throw/catch
        self.basedir   = basedir

        self.runname   = runsetname
        self.rundir    = os.path.join(basedir,runsetname)
        if os.path.isdir(self.rundir):
            print "Oops! Already such-named run in this base folder!"
            sys.exit(0) #TODO: Should prob. use some throw/catch
        os.mkdir(self.rundir);

    def newRun(self,method,runname,N):
        """
        Setup and do a new run

        Inputs:
         - method:  Which method do we use? (see run.__init__ for possible values)
         - runname: Name of the run
         - N:       Number of points

         Saves the created run object to the runs dictionary (=HashMap in Java),
         and also returns it.
        """
        thisrun = run(self.rundir, runname,
                      self.m,self.l,self.omega,self.A,
                      self.viscosity,self.phi_0,self.v_0,
                      self.t_end, N, method,self.constructRunlabel)

        thisrun.doRun()
        
        self.runs[self.constructRunlabel(method,runname)] = thisrun
        return thisrun

    def getRun(self,method,runname):
        return self.runs[constructRunlabel(method,runname)]
    def getRun(self,dataset):
        return self.runs[dataset]
    
    def constructRunlabel(self,method,runname):
        """
        Returns a string containing the filename/hashmap-index
        corresponding to method and runname
        """
        return (method + "-" + str(runname))
        
class run:
    """
    This class contains the info for a single run.

    The outdata is stored as a file with format \"time pos dpos/dt\".
    Path of the file is (on UNIX) /folder/method-runname.out

    Aviable methods (from pendulum):
     euler, euler_cromer, midpoint, euler_richardson, half_step, rk2, rk4, asc
    """
    def __init__(self,folder,runname,m,l,omega,A,viscosity,phi_0,v_0,t_end, N, method,runlabel_constructor):
        """
        Constructor: setup class variables

        Inputs (the non-obvious ones):

         - folder:  Where to store outdata
         - runname: Name of this run
         - method:  Which solver to run
         - t_end:   Final time AS A MULTIPLUM OF PI
         - runlabel_constructor: Method that constructs rulabels
        """
        self.folder    = str(folder)
        self.runname   = str(runname)
        self.m         = float(m)
        self.l         = float(l)
        self.omega     = float(omega)
        self.A         = float(A)
        self.viscosity = float(viscosity)
        self.phi_0     = float(phi_0)
        self.v_0       = float(v_0)
        self.t_end     = float(t_end)
        self.N         = int(N)
        self.method    = str(method)

        self.runfile   = os.path.join(self.folder,runlabel_constructor(method,runname) + ".out")

    def doRun(self):
        """Run a simulation"""
        pendObj = pendelum(self.m,self.l,self.omega, self.A,
                           self.viscosity, self.phi_0, self.v_0,
                           self.N, self.t_end)
        if   self.method == "euler":
            pendObj.euler();
        elif self.method == "euler_cromer":
            pendObj.euler_cromer()
        elif self.method == "midpoint":
            pendObj.midpoint()
        elif self.method == "euler_richardson":
            pendObj.euler_richardson()
        elif self.method == "half_step":
            pendObj.half_step()
        elif self.method == "rk2":
            pendObj.rk2()
        elif self.method == "rk4":
            pendObj.rk4()
        elif self.method == "asc":
            pendObj.asc() #This method doesn't seem to work!
        else:
            print "Invalid method!";
            system.exit(0) #TODO: Should prob. use some throw/catch

        #Move the output file to the correct spot
        shutil.copyfile(self.method + ".out", self.runfile)
        os.remove(self.method + ".out")

    def plotTheta(self):
        """Use plotFile to plot Theta(t)"""
        plot = plotFile(self.runfile,"","","1:2")
        del plot
    def plotdThetadt(self):
        """Use plotFile to plot dTheta/dt(t)"""
        plot = plotFile(self.runfile,"","","1:3")
        del plot
        
class plotFile:
    """
    This is used to setup gnuplot and plot the contents of one file
    """

    def __init__(self,plotfile,xlabel,ylabel,using):
        """
        Setup a plot (create a gnuplot script)

        Input:
         - plotfile: File to plot
         - xlabel:   Label of x-axis
         - ylabel:   Label of y-axis
         - using:    Which columns do we plot (ex. \"1:2\" to plot
                     column 2 vs column 1)

        """

        self.ofilename = plotfile + ".gnuplot"

        ofile = open(self.ofilename,'w')
        ofile.write("plot '" + plotfile + "' using " + using + " w l\n")
        ofile.close()
        cmd = "gnuplot -persist " + self.ofilename
        import commands
        commands.getstatusoutput(cmd)
        
    def __del__(self):
        """
        Destructor: Deletes the plotfile
        """
        os.remove(self.ofilename)
