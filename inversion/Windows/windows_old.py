import numpy as np
import os
import shutil
from scipy import signal
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess


class windows:
    def __init__(self,par):
        self.par = par
        # PARMAMETERS

        # Reading Windows Directory 
        self.path = self.par.InvPar['-WinPath'] + 'WIN/'
        os.chdir(self.path)
        command = "ls -l *N.win.qual | awk '{print $9}' | awk -F- '{print $1}'"
        self.stat_name = subprocess.check_output(command,shell=True).split('\n')[:-1]

        # Windows number for every trace
        self.nw = np.zeros((len(self.stat_name),3))

        # Number of windows
        self.stat_win = {}

        # Windows Information [timeI,timeF,Tshift,CC,dlnA]
        self.stat_we = {}
        self.stat_wn = {}
        self.stat_wz = {}


        # Statistical Array

        self._win_e = np.zeros((1,5))
        self._win_n = np.zeros((1,5))
        self._win_z = np.zeros((1,5))

        
        self.CreateWindows()

        # Reading WorkDir
        os.chdir(self.par.InvPar["-WorkDir"])

    def CreateWindows(self):

        i = 0
        for istat in self.stat_name:
            win_e = istat + "-HHE.win.qual"
            win_n = istat + "-HHN.win.qual"
            win_z = istat + "-HHZ.win.qual"
            
            win_all = np.array([nWin(win_e),nWin(win_n),nWin(win_z)])
            self.stat_win.update({istat:win_all})

            self.Wins(istat)

        
            
    def Wins(self,istat):
        win_e = istat + "-HHE.win.qual"
        win_n = istat + "-HHN.win.qual"
        win_z = istat + "-HHZ.win.qual"
            
        # EAST COMPONENT
        if self.stat_win[istat][0] > 0:
            inf = ReadWin(win_e,self.stat_win[istat][0])
            self._win_e = np.append(self._win_e,inf,axis=0)
            self.stat_we.update({istat:inf})

        # NORTH COMPONENT
        if self.stat_win[istat][1] > 0:
            inf = ReadWin(win_n,self.stat_win[istat][1])
            self._win_n = np.append(self._win_n,inf,axis=0)
            self.stat_wn.update({istat:inf})

        # EAST COMPONENT
        if self.stat_win[istat][2] > 0:
            inf = ReadWin(win_z,self.stat_win[istat][2])
            self._win_z = np.append(self._win_z,inf,axis=0)         
            self.stat_wz.update({istat:inf})


    def Hist(self):


        # Delete First Empty Row

        _win_e = np.delete(self._win_e,0,0)
        _win_n = np.delete(self._win_n,0,0)
        _win_z = np.delete(self._win_z,0,0)

        
        plt.figure()
        plt.subplot(1,3,1)
        plt.title(r'$\tau\ [s]$')
        plt.hist(_win_e[:,2],color='k')
        plt.subplot(1,3,2)
        plt.title(r'$CC$')
        plt.hist(_win_n[:,3],color='k')
        plt.subplot(1,3,3)
        plt.title(r'$\Delta\ lnA$')
        plt.hist(_win_z[:,4],color='k')
        plt.show()
        

        
        
            
def nWin(name):
    f = open(name,'r')
    nwin = int(f.readline().split('=')[1])
    f.close()
    return nwin

def ReadWin(name,nwin):
    f = open(name,'r')

    # Ignore Header
    f.readline()
    f.readline()

    i = 0
    t_win = np.zeros([nwin,5])
    while i < nwin:
        line = f.readline().strip(' ').strip('\n').split()
        t_win[i,:] = np.array([float(line[1]),float(line[2]),float(line[3]), \
                               float(line[4]),float(line[5])])
        
        i +=1

    return t_win

    f.close()
