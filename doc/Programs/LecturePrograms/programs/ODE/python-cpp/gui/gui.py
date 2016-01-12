#!/usr/bin/env python
import runcontrol

import Tkinter,Pmw

class control:
    runset = None;

    def __init__(self):
        wd_main = Tkinter.Tk()
        wd_main.title("Pendulum RunControl GUI")
        self.wd_main_frame = Tkinter.Frame(wd_main)
        self.wd_main_frame.pack()

        configpane(self.wd_main_frame,self)

        wd_main.mainloop()
        
    def createRunset(self,folder,name,m,l,omega,A,viscosity,phi0,v0,tend):
        self.runset = runcontrol.runset(folder,name,m,l,omega,A,viscosity,phi0,v0,tend)
        runpane(self.wd_main_frame,self)
        
    def run(self,method,runname,N):
        self.runset.newRun(method,runname,N)

    def constructRunlabel(self,method,runname):
        """
        Returns a string containing the filename/hashmap-index
        corresponding to method and runname
        """
        return self.runset.constructRunlabel(method,runname)
    
    def plotFunc(self,dataset,what):
        if what == "theta":
            self.runset.getRun(dataset).plotTheta()
        elif what == "dtheta/dt":
            self.runset.getRun(dataset).plotdThetadt()
        else:
            print "WHAT?"
            system.exit(0)

    def compare(self,dataset1, dataset2):
        """
        Makes a comparison of theta(t) for plot1 and plot2

        Input:
        - dataset1, dataset2: Names of the datasets

        No interpolation: N must be equal for both datasets!
        """
        #Get the filenames and files
        run1 = self.runset.getRun(dataset1)
        run2 = self.runset.getRun(dataset2)
        runfile1 = self.runset.getRun(dataset1).runfile
        runfile2 = self.runset.getRun(dataset2).runfile

        rf1 = open(runfile1,'r')
        rf2 = open(runfile2,'r')

        #Create a filename for writing and open it
        comparename = "compare-"+dataset1+"-VS-"+dataset2
        import os
        outname = os.path.join(self.runset.rundir,comparename)
        of = open(outname,'w')

        #Create an outfile in format "time theta1(t)-theta2(t) log(abs(deltatheta(t)))
        import math
        while True:
            lrf1 = rf1.readline()
            lrf2 = rf2.readline()

            if not lrf1 or not lrf2:
                break;

            [str_t1,str_theta1,str_v1] = lrf1.split()
            [str_t2,str_theta2,str_v2] = lrf2.split()
            t1     = float(str_t1)
            t2     = float(str_t2)
            theta1 = float(str_theta1)
            theta2 = float(str_theta2)
            v1     = float(str_v1)
            v2     = float(str_v2)
            #Check for desync
            if t1 != t2:
                print "Desync!"
                return
            #Calculate
            dt = theta1-theta2
            logdt = math.log(math.fabs(theta1-theta2))

            #Write to file
            of.write("%30.20E %30.20E %30.20E\n" % (t1,dt,logdt))

        #Close files
        rf1.close()
        rf2.close()
        of.close()

        #Plot dt(t):
        plot1 = runcontrol.plotFile(outname,"","","1:2")
        del plot1
        plot2 = runcontrol.plotFile(outname,"","","1:3")
        del plot2
        
        
class configpane:
    
    disablelist = []
    controller  = None #Used for callback
    
    def __init__(self,parent,controller):

        self.controller = controller

        #Create configpane
        parentframe = Tkinter.Frame(parent)
        parentframe.pack(side='left')
        
        parentframe_label = Tkinter.Label(parentframe,text="Configure runset:")
        parentframe_label.grid(row=0, column=0, columnspan=2)
        
        parentframe_folderlabel = Tkinter.Label(parentframe,text="Save dataset in folder:")
        parentframe_folderlabel.grid(row=1,column=0)
        self.str_folder = Tkinter.StringVar()
        self.str_folder.set("/tmp")
        parentframe_folderentry = Tkinter.Entry(parentframe,textvariable=self.str_folder)
        parentframe_folderentry.grid(row=1,column=1)
        self.disablelist.append(parentframe_folderentry)
        
        parentframe_mlabel = Tkinter.Label(parentframe,text="Mass m:")
        parentframe_mlabel.grid(row=2,column=0)
        self.str_m = Tkinter.StringVar()
        self.str_m.set("1.0")
        parentframe_mentry = Tkinter.Entry(parentframe,textvariable=self.str_m)
        parentframe_mentry.grid(row=2,column=1)
        self.disablelist.append(parentframe_mentry)
                
        parentframe_label = Tkinter.Label(parentframe,text="Length l:")
        parentframe_label.grid(row=3,column=0)
        self.str_l = Tkinter.StringVar()
        self.str_l.set("1.0")
        parentframe_lentry = Tkinter.Entry(parentframe,textvariable=self.str_l)
        parentframe_lentry.grid(row=3,column=1)
        self.disablelist.append(parentframe_lentry)
        
        parentframe_omegalabel = Tkinter.Label(parentframe,text="Omega of the force:")
        parentframe_omegalabel.grid(row=4,column=0)
        self.str_omega = Tkinter.StringVar()
        self.str_omega.set("1.0")
        parentframe_omegaentry = Tkinter.Entry(parentframe,textvariable=self.str_omega)
        parentframe_omegaentry.grid(row=4,column=1)
        self.disablelist.append(parentframe_omegaentry)        
        
        parentframe_Alabel = Tkinter.Label(parentframe,text="Amplitude of the force:")
        parentframe_Alabel.grid(row=5,column=0)
        self.str_A = Tkinter.StringVar()
        self.str_A.set("1.0")
        parentframe_Aentry = Tkinter.Entry(parentframe,textvariable=self.str_A)
        parentframe_Aentry.grid(row=5,column=1)
        self.disablelist.append(parentframe_Aentry)
        
        parentframe_viscositylabel = Tkinter.Label(parentframe,text="Viscous drag value (viscosity):")
        parentframe_viscositylabel.grid(row=6,column=0)
        self.str_viscosity = Tkinter.StringVar()
        self.str_viscosity.set("1.0")
        parentframe_viscosityentry = Tkinter.Entry(parentframe,textvariable=self.str_viscosity)
        parentframe_viscosityentry.grid(row=6,column=1)
        self.disablelist.append(parentframe_viscosityentry)
        
        parentframe_phi0label = Tkinter.Label(parentframe,text="Phi0:")
        parentframe_phi0label.grid(row=7,column=0)
        self.str_phi0 = Tkinter.StringVar()
        self.str_phi0.set("2.0")
        parentframe_phi0entry = Tkinter.Entry(parentframe,textvariable=self.str_phi0)
        parentframe_phi0entry.grid(row=7,column=1)
        self.disablelist.append(parentframe_phi0entry)
        
        parentframe_v0label = Tkinter.Label(parentframe,text="v0:")
        parentframe_v0label.grid(row=8,column=0)
        self.str_v0 = Tkinter.StringVar()
        self.str_v0.set("0.0")
        parentframe_v0entry = Tkinter.Entry(parentframe,textvariable=self.str_v0)
        parentframe_v0entry.grid(row=8,column=1)
        self.disablelist.append(parentframe_v0entry)
        
        parentframe_tendlabel = Tkinter.Label(parentframe,text="Final time step (as multiplum of PI):")
        parentframe_tendlabel.grid(row=9,column=0)
        self.str_tend = Tkinter.StringVar()
        self.str_tend.set("4.0")
        parentframe_tendentry = Tkinter.Entry(parentframe,textvariable=self.str_tend)
        parentframe_tendentry.grid(row=9,column=1)
        self.disablelist.append(parentframe_tendentry)
        
        parentframe_namelabel = Tkinter.Label(parentframe,text="Name of runset:")
        parentframe_namelabel.grid(row=10,column=0)
        self.str_name = Tkinter.StringVar()
        self.str_name.set("testset")
        parentframe_nameentry = Tkinter.Entry(parentframe,textvariable=self.str_name)
        parentframe_nameentry.grid(row=10,column=1)
        self.disablelist.append(parentframe_nameentry)

        parentframe_createbutton = Tkinter.Button(parentframe, text='Create new runset',command=self.createRunset)
        parentframe_createbutton.grid(row=11,column=0,columnspan=2)
        self.disablelist.append(parentframe_createbutton)
        
    def createRunset(self):
        """
        Call the controller with the correct arguments and then lockdown
        """
        import os
        path = os.path.join(self.str_folder.get(),self.str_name.get())
        if os.path.isdir(path):
            import tkMessageBox
            if tkMessageBox.askokcancel("Already exists","A run with this name is already saved here. Overwrite?"):
                import shutil
                shutil.rmtree(path)
                tkMessageBox.Message(type="ok",message=("Folder " + path + " was removed")).show()
            else:
                return

        self.controller.createRunset(self.str_folder.get(),self.str_name.get(),
                                     float(self.str_m.get()),float(self.str_l.get()), float(self.str_omega.get()),
                                     float(self.str_A.get()),float(self.str_viscosity.get()),
                                     float(self.str_phi0.get()),float(self.str_v0.get()),float(self.str_tend.get()))
        self.lockdown()

    def lockdown(self):
        """
        Grays out all text entry boxes, and the button
        """
        for widget in self.disablelist:
            widget.config(state='disabled')

class runpane:

    #Object variables initialization
    compare_plot1 = None;
    
    def __init__(self,parent,controller):
        self.controller = controller
        
        parentframe = Tkinter.Frame(parent)
        parentframe.pack(side='left')

        parentframe_newrunlabel = Tkinter.Label(parentframe,text="New run:")
        parentframe_newrunlabel.grid(row=0,column=0,columnspan=2)

        parentframe_methodlabel = Tkinter.Label(parentframe,text="Method:")
        parentframe_methodlabel.grid(row=1,column=0)
        self.str_method = Tkinter.StringVar()
        self.str_method.set("euler")
        parentframe_methoddropdown = Pmw.OptionMenu(parentframe,
                                                    items=['euler','euler_cromer','midpoint','euler_richardson','half_step','rk2','rk4'],
                                                    menubutton_textvariable = self.str_method,
                                                    menubutton_width = 10)
        parentframe_methoddropdown.grid(row=1,column=1)

        parentframe_Nlabel = Tkinter.Label(parentframe,text="N:")
        parentframe_Nlabel.grid(row=2,column=0)
        self.str_N = Tkinter.StringVar()
        self.str_N.set("10000")
        parentframe_Nentry = Tkinter.Entry(parentframe,textvariable=self.str_N)
        parentframe_Nentry.grid(row=2,column=1)

        parentframe_namelabel = Tkinter.Label(parentframe,text="Runname:")
        parentframe_namelabel.grid(row=3,column=0)
        self.str_name = Tkinter.StringVar()
        self.update_str_name()
        #parentframe_Nentry.bind('<Return>',self.update_str_name)
        parentframe_nameentry = Tkinter.Entry(parentframe,textvariable=self.str_name,)
        parentframe_nameentry.grid(row=3,column=1)
        parentframe_nameupdatebutton = Tkinter.Button(parentframe, text="Update runname",command=self.update_str_name)
        parentframe_nameupdatebutton.grid(row=4,column=0,columnspan=2)

        parentframe_runbutton = Tkinter.Button(parentframe, text = "RUN!", command=self.run)
        parentframe_runbutton.grid(row=5,column=0,columnspan = 2)

        self.parentframe_conductedrunslist = Pmw.ScrolledListBox(parentframe,
                                                                 listbox_selectmode="single",
                                                                 vscrollmode="static", hscrollmode="none",
                                                                 label_text="Conducted runs:",labelpos='n',
                                                                 listbox_width=15,listbox_height=10)
        self.parentframe_conductedrunslist.grid(row=6,column=0,columnspan=2)

        parentframe_plot1button = Tkinter.Button(parentframe, text="Plot theta(t)",command=self.plotTheta)
        parentframe_plot1button.grid(row=7,column=0)

        parentframe_plot2button = Tkinter.Button(parentframe, text="Plot dtheta/dt(t)",command=self.plotV)
        parentframe_plot2button.grid(row=7,column=1)

        parentframe_comparebutton = Tkinter.Button(parentframe, text="Compare theta(t)...",command=self.compare)
        parentframe_comparebutton.grid(row=8,column=0,columnspan=2)
        
    def update_str_name(self,event=None):
        self.str_name.set("N="+self.str_N.get())

    def run(self):
        try:
            N = int(self.str_N.get())
        except ValueError:
            import tkMessageBox
            tkMessageBox.Message(type="ok",message=("Please enter a number in the \"N\" box")).show()
            return

        #Search through the listbox and check that it doesnt already exist
        list = self.parentframe_conductedrunslist.get()
        for item in list:
            if item == self.controller.constructRunlabel(self.str_method.get(),self.str_name.get()):
                import tkMessageBox
                tkMessageBox.Message(type="ok",message=("This run already exists, please select another name")).show()
                return
        
        self.controller.run(self.str_method.get(),self.str_name.get(),N)

        self.parentframe_conductedrunslist.insert('end', self.controller.constructRunlabel(self.str_method.get(),self.str_name.get()))
    def plotTheta (self):
        try:
            index = self.parentframe_conductedrunslist.getvalue()[0];
        except IndexError:
            import tkMessageBox
            tkMessageBox.Message(type="ok",message=("Please select a run to plot")).show()
            return
        self.controller.plotFunc(index,"theta")
    
    def plotV(self):
        try:
            index = self.parentframe_conductedrunslist.getvalue()[0];
        except IndexError:
            import tkMessageBox
            tkMessageBox.Message(type="ok",message=("Please select a run to plot")).show()
            return
        self.controller.plotFunc(index,"dtheta/dt")

    def compare(self):
        """Setup a comparison of two plots"""
        try:
            index = self.parentframe_conductedrunslist.getvalue()[0];
        except IndexError:
            import tkMessageBox
            tkMessageBox.Message(type="ok",message=("Please select a run to compare")).show()
            return

        if self.compare_plot1 == None:
            #First push
            self.compare_plot1 = index
            import tkMessageBox
            tkMessageBox.Message(type="ok",message=("Select another run to compare with, and push \"Compare...\" again")).show()
        else:
            #Secound push
            compare_plot2 = index
            self.controller.compare(self.compare_plot1,compare_plot2)
            self.compare_plot1 = None
        
        
        

### MAIN PROGRAM ###

#Setup main window with dataentry pane


controller = control()

