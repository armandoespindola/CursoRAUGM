import numpy as np
import os
from inversion.Tools.UtilFunc import butter_bandpass_filter,PlotNorm,butter_lowpass_filter
import math
#from inversion.Windows.windows import windows
from scipy import signal

class ForwardSource:
    def __init__(self,par):
        self.par = par
        self.filter = butter_bandpass_filter
        self.fw = np.fft.rfftfreq(self.par.nt,self.par.dt)
        
    def ReadDir(self,i):
        os.chdir(self.par.InvPar["-WorkDir"])
        os.chdir(self.par.EventName[i] + "/DATA")
        

    def FwdSrc(self):

        freq = self.par.freq[-1]
        
        for k in range(0,self.par.EventNum):
            os.chdir(self.par.InvPar["-WorkDir"])
            os.chdir(self.par.EventName[k])

            #src = self.filt(self.filt(np.fromfile('SOURCE.bin',dtype='f4'),freq)[::-1],freq)[::-1]
            src = np.fromfile('SOURCE.bin',dtype='f4')
            src = self.filt(np.fromfile('SOURCE.bin',dtype='f4'),self.par.freq)
            src_filt = src
            #amp_s,ph_s = FreqExtract(src)
            #f0 = 1.0 / (self.par.dt * self.par.nt)
            #fact = 2.0 / self.par.nt
            #for i in range(0,self.par.freq.size):
            #    nf = int(round(self.par.freq[i] / f0))
            #    for nt in range(0,self.par.nt-1):
            #        arg = 2.0 * np.pi * self.fw[nf] * self.par.dt * nt
            #        pcosx = np.cos(arg + ph_s[nf])
            #        src_filt[nt] = src_filt[nt] + amp_s[nf] * (pcosx) * fact

            src_filt = Taper(src_filt)

            #src_filt = src #FreqExtract(src,self.par.nt,self.par.dt,self.par.freq)
            src_filt.astype('f4').tofile('SrcTime.bin')

    def filt(self,data,freq):
        data_filt = self.filter(data,freq[0],freq[-1],1.0/self.par.dt,order=4)
        return data_filt



def Taper(data):
    win = signal.tukey(data.size,alpha=0.2)
    return win * data



def FreqExtract(data):
    dataf = np.fft.rfft(data)
    phase = np.unwrap(np.angle(dataf))
    amp = np.abs(dataf)
    return [amp,phase]
