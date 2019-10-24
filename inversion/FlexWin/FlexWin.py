import numpy as np
import os
import shutil
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from scipy import signal
from obspy import read,io
from Tools.UtilFunc import butter_bandpass_filter,PlotNorm,butter_lowpass_filter

class FlexWin:
    def __init__(self,par):
        # PARAMETERS
        self.par = par
        self.filter = butter_bandpass_filter

    def Compile(self):
	os.chdir(self.par.InvPar["-WorkDir"])
	os.chdir("flexwin/")
        subprocess.call("make -f make_gfortran clean ",shell=True)
        subprocess.call("make -f make_gfortran ",shell=True)

    def CopySolver(self):
        os.chdir(self.par.InvPar["-WorkDir"])
	os.chdir("flexwin/")
	for i in range(0,self.par.EventNum):
	    dest = self.par.InvPar["-WorkDir"]+"/"+self.par.EventName[i]+ \
	           "/"+self.par.EventName[i]+"-FlexWin.out"
	    shutil.copy("flexwin",dest)

    def CreateInputFlexWin(self):
        # Period For Windows
        T2 = 1.0 / self.par.freq[0]
        T1 = 1.0 / self.par.freq[-1]
        print self.par.freq[0],self.par.freq[-1]
        for i in range(0,self.par.EventNum):
            os.chdir(self.par.InvPar["-WorkDir"])
            os.chdir(self.par.EventName[i])
            with open('PAR_FILE', 'r') as file:
                data = file.readlines()
                data[23] = "WIN_MIN_PERIOD                  = " + str(T1) +"\n"
                data[24] = "WIN_MAX_PERIOD                  = " + str(T2) + "\n"
            with open('PAR_FILE', 'w') as file:
                file.writelines( data )

            # Writing Input File FlexWin
            command = "ls -l DATA/*VX.bin | awk '{print $9}' | awk -F- '{print $1}' | awk -F/ '{print $2}'"
            stat =subprocess.check_output(command,shell=True).split('\n')
            f = open("FlexWin.input","w")
            ntrace = (len(stat)-1) * 3
            f.write('%d\n' % (ntrace))
            for i in range(0,len(stat)-1):
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHE.OBS.SAC\n")
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHE.SYN.SAC\n")
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHE\n")
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHN.OBS.SAC\n")
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHN.SYN.SAC\n")
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHN\n")
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHZ.OBS.SAC\n")
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHZ.SYN.SAC\n")
                f.write(self.par.InvPar["-WinPath"] + "WIN/" + stat[i] + "-HHZ\n")
            f.close()





    def CreateSacFiles(self,veldisp,freqmin):
        for i in range(0,self.par.EventNum):
            os.chdir(self.par.InvPar["-WorkDir"])
            os.chdir(self.par.EventName[i])
            #Station
            stat = np.genfromtxt('receptors.par',dtype='str',comments='#',delimiter=',')
            n = stat.shape # n[0] number of stations that record the event
            subprocess.call("mkdir -p " + self.par.InvPar["-WinPath"] + "WIN",shell=True)
            os.chdir('DATA/')
            for i in range(0,n[0]):
                # Syntethic Data
                sx_v = np.fromfile(stat[i,0].strip()+'-VX.bin',dtype='f4')
                sy_v = np.fromfile(stat[i,0].strip()+'-VY.bin',dtype='f4')
                sz_v = np.fromfile(stat[i,0].strip()+'-VZ.bin',dtype='f4')
                # Observed Data
                o_x_v = read(stat[i,0] + "-HHE-OBS.sac")
                o_y_v = read(stat[i,0] + "-HHN-OBS.sac")
                o_z_v = read(stat[i,0] + "-HHZ-OBS.sac")

                o_x =  o_x_v.copy()
                o_x[0].data[:] = 0.0
                o_y =  o_y_v.copy()
                o_y[0].data[:] = 0.0
                o_z =  o_z_v.copy()
                o_z[0].data[:] = 0.0



                if veldisp == 'disp' :
                    sx = np.cumsum(sx_v) * self.par.dt
                    sy = np.cumsum(sy_v) * self.par.dt
                    sz = np.cumsum(sz_v) * self.par.dt

                    o_x[0].data = np.cumsum(o_x_v.data[:]) * self.par.dt
                    o_y[0].data = np.cumsum(o_y_v.data[:]) * self.par.dt
                    o_z[0].data = np.cumsum(o_z_v.data[:]) * self.par.dt

                elif veldisp == 'vel':
                    sx = sx_v
                    sy = sy_v
                    sz = sz_v

                    o_x[0].data = o_x_v[0].data
                    o_y[0].data = o_y_v[0].data
                    o_z[0].data = o_z_v[0].data
                    

                # Filter data and Create Input File Flexwin
                s_x = o_x.copy()
                s_x[0].data[:] = 0.0

                # s_x[0].data[:len(sx)] = self.filter(self.filter(Taper(sx[::1],0.05),self.par.freq[-1] \
                #                                                 ,1.0/self.par.dt, order=2)[::-1], \
                #                                     self.par.freq[-1],1.0/self.par.dt, order=2)[::-1]
                
                # o_x[0].data = self.filter(self.filter(Taper(o_x[0].data,0.05),self.par.freq[-1] \
                #                                       ,1.0/self.par.dt, order=2)[::-1], \
                #                           self.par.freq[-1],1.0/self.par.dt, order=2)[::-1]
                

                s_x[0].data[:len(sx)] = self.filter(Taper(sx[::1],0.05),self.par.freq[0],self.par.freq[-1] \
                                                                ,1.0/self.par.dt, order=2)
                
                o_x[0].data = self.filter(Taper(o_x[0].data,0.05),self.par.freq[0],self.par.freq[-1] \
                                                      ,1.0/self.par.dt, order=2)


                s_y = o_y.copy()
                s_y[0].data[:] = 0.0
                
                # s_y[0].data[:len(sy)] = self.filter(self.filter(Taper(sy[::1],0.05),self.par.freq[-1] \
                #                                                 ,1.0/self.par.dt, order=2)[::-1], \
                #                                     self.par.freq[-1],1.0/self.par.dt, order=2)[::-1]
                
                # o_y[0].data = self.filter(self.filter(Taper(o_y[0].data,0.05),self.par.freq[-1] \
                #                                       ,1.0/self.par.dt, order=2)[::-1], \
                #                           self.par.freq[-1],1.0/self.par.dt, order=2)[::-1]


                s_y[0].data[:len(sy)] = self.filter(Taper(sy[::1],0.05),self.par.freq[0],self.par.freq[-1] \
                                                                ,1.0/self.par.dt, order=2)
                
                o_y[0].data = self.filter(Taper(o_y[0].data,0.05),self.par.freq[0],self.par.freq[-1] \
                                                      ,1.0/self.par.dt, order=2)
                
                s_z = o_z.copy()
                s_z[0].data[:] = 0.0

                # s_z[0].data[:len(sz)] = self.filter(self.filter(Taper(sz[::1],0.05),self.par.freq[-1] \
                #                                                 ,1.0/self.par.dt, order=2)[::-1], \
                #                                     self.par.freq[-1],1.0/self.par.dt, order=2)[::-1]

                # o_z[0].data = self.filter(self.filter(Taper(o_z[0].data,0.05),self.par.freq[-1] \
                #                                       ,1.0/self.par.dt, order=2)[::-1], \
                #                           self.par.freq[-1],1.0/self.par.dt, order=2)[::-1]

                s_z[0].data[:len(sz)] = self.filter(Taper(sz[::1],0.05),self.par.freq[0],self.par.freq[-1] \
                                                                ,1.0/self.par.dt, order=2)

                o_z[0].data = self.filter(Taper(o_z[0].data,0.05),self.par.freq[0],self.par.freq[-1] \
                                                      ,1.0/self.par.dt, order=2)



                o_x[0].write("WIN/" + stat[i,0] + "-HHE.OBS.SAC",format="SAC")
                s_x[0].write("WIN/" + stat[i,0] + "-HHE.SYN.SAC",format="SAC")

                o_y[0].write("WIN/" + stat[i,0] + "-HHN.OBS.SAC",format="SAC")
                s_y[0].write("WIN/" + stat[i,0] + "-HHN.SYN.SAC",format="SAC")

                o_z[0].write("WIN/" + stat[i,0] + "-HHZ.OBS.SAC",format="SAC")
                s_z[0].write("WIN/" + stat[i,0] + "-HHZ.SYN.SAC",format="SAC")
                
        
    def RunFlexWin(self):
	os.chdir(self.par.InvPar["-WorkDir"]);
	for i in range(0,self.par.EventNum):
	    program = "mpirun -np " + "1"+ \
		      " ./" + self.par.EventName[i] +"-FlexWin.out" \
		      + " < " +self.par.InvPar["-FlexWinInput"] + \
		      " > FlexWin.log";
	    print "Running (FlexWin): " + self.par.EventName[i] +" -FlexWin.out"
	    os.chdir(self.par.EventName[i]);
            subprocess.call("ln -sf ../flexwin/iasp91.hed",shell=True);
            subprocess.call("ln -sf ../flexwin/iasp91.tbl",shell=True);
	    subprocess.call(program,shell=True);
	    os.chdir("../")



def Taper(data,alpha=0.25):
    win = signal.tukey(data.size,alpha)
    return win * data
