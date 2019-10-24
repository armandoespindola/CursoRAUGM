import numpy as np
import os
from inversion.Tools.UtilFunc import butter_bandpass_filter,PlotNorm,butter_lowpass_filter
import math
from inversion.Windows.windows import windows
from scipy import signal

class AdjointSource:
    def __init__(self,par):
        self.par = par
        print '## freq: ',self.par.freq
        self.filter = butter_bandpass_filter
        
    def ReadDir(self,i):
        os.chdir(self.par.InvPar["-WorkDir"])
        os.chdir(self.par.EventName[i] + "/DATA")
        

    def AdjSrc(self,ntype,veldisp,iter):

        #freq = self.par.freq[-1] #+ (2.0 / self.par.Tsim)
        freq = self.par.freq[-1]
        print freq
        
        for k in range(0,self.par.EventNum):
            os.chdir(self.par.InvPar["-WorkDir"])
            os.chdir(self.par.EventName[k])
            # Call Windows Class
            win = windows(self.par)
            
            os.chdir(self.par.EventName[k])
            print self.par.EventName[k]

            #Station
            stat = np.genfromtxt('receptors.par',dtype='str',comments='#',delimiter=',')
            n = stat.shape # n[0] number of stations that record the event

            os.chdir('DATA/')
            
            for i in range(0,n[0]):
                rx_v = np.fromfile(stat[i,0].strip()+'-VX.bin',dtype='f4')
                ry_v = np.fromfile(stat[i,0].strip()+'-VY.bin',dtype='f4')
                rz_v = np.fromfile(stat[i,0].strip()+'-VZ.bin',dtype='f4')
                obs_rx_v = np.fromfile(stat[i,0].strip()+'-VX-OBS.bin',dtype='f4')
                obs_ry_v = np.fromfile(stat[i,0].strip()+'-VY-OBS.bin',dtype='f4')
                obs_rz_v = np.fromfile(stat[i,0].strip()+'-VZ-OBS.bin',dtype='f4')

                if veldisp == 'disp' :
                    rx = signal.detrend(np.cumsum(rx_v) * self.par.dt)
                    ry = signal.detrend(np.cumsum(ry_v) * self.par.dt)
                    rz = signal.detrend(np.cumsum(rz_v) * self.par.dt)

                    obs_rx = signal.detrend(np.cumsum(obs_rx_v) * self.par.dt)
                    obs_ry = signal.detrend(np.cumsum(obs_ry_v) * self.par.dt)
                    obs_rz = signal.detrend(np.cumsum(obs_rz_v) * self.par.dt)

                elif veldisp == 'vel':
                    rx = rx_v
                    ry = ry_v
                    rz = rz_v

                    obs_rx = obs_rx_v
                    obs_ry = obs_ry_v
                    obs_rz = obs_rz_v



                # Filter Seismograms - Zero Phase

                rx = self.filt(Taper(rx,0.1),self.par.freq)
                ry = self.filt(Taper(ry,0.1),self.par.freq)
                rz = self.filt(Taper(rz,0.1),self.par.freq)

                obs_rx = self.filt(Taper(obs_rx,0.1),self.par.freq)
                obs_ry = self.filt(Taper(obs_ry,0.1),self.par.freq)
                obs_rz = self.filt(Taper(obs_rz,0.1),self.par.freq)


                # Initialize Adjoint Variables
                adj_rx = rx * 0.0
                adj_ry = ry * 0.0
                adj_rz = rz * 0.0
                
                stat_name = stat[i,0]
                # Number of windows [east,north,z]
                nwin = win.stat_win[stat_name]

                # East Component 
                if nwin[0] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    inf_w = win.stat_we[stat_name]
                    for i_nw in range(0,nwin[0]):
                        nt_start = int(np.floor(inf_w[i_nw,0] / self.par.dt))
                        nt_end = int(np.floor(inf_w[i_nw,1] / self.par.dt))
                        d_syn = Taper(rx[nt_start:nt_end],0.05)
                        d_obs = Taper(obs_rx[nt_start:nt_end],0.05)
                        adj_rx[nt_start:nt_end] = self.Adj(d_syn,d_obs,ntype)
                    adj_rx = self.filt(self.filt(adj_rx,self.par.freq,order=1)[::-1],self.par.freq,order=1)[::-1]
                    adj_rx.astype('f4').tofile(stat[i,0].strip()+'-ADJX.bin')
                elif nwin[0] == 0:
                    adj_rx.astype('f4').tofile(stat[i,0].strip()+'-ADJX.bin')
                elif nwin[0] < 0:
                    adj_rx = self.Adj(rx,obs_rx,ntype)
                    adj_rx.astype('f4').tofile(stat[i,0].strip()+'-ADJX.bin')



                # North Component 
                if nwin[1] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    inf_w = win.stat_wn[stat_name]
                    for i_nw in range(0,nwin[1]):
                        nt_start = int(np.floor(inf_w[i_nw,0] / self.par.dt))
                        nt_end = int(np.floor(inf_w[i_nw,1] / self.par.dt))
                        d_syn = Taper(ry[nt_start:nt_end],0.05)
                        d_obs = Taper(obs_ry[nt_start:nt_end],0.05)
                        adj_ry[nt_start:nt_end] = self.Adj(d_syn,d_obs,ntype)
                    adj_ry = self.filt(self.filt(adj_ry,self.par.freq,order=1)[::-1],self.par.freq,order=1)[::-1]
                    adj_ry.astype('f4').tofile(stat[i,0].strip()+'-ADJY.bin')
                elif nwin[1] == 0:
                    adj_ry.astype('f4').tofile(stat[i,0].strip()+'-ADJY.bin')
                elif nwin[1] < 0:
                    adj_ry = self.Adj(ry,obs_ry,ntype)
                    adj_ry.astype('f4').tofile(stat[i,0].strip()+'-ADJY.bin')


                # Z Component 
                if nwin[2] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    inf_w = win.stat_wz[stat_name]
                    for i_nw in range(0,nwin[2]):
                        nt_start = int(np.floor(inf_w[i_nw,0] / self.par.dt))
                        nt_end = int(np.floor(inf_w[i_nw,1] / self.par.dt))
                        d_syn = Taper(rz[nt_start:nt_end],0.05)
                        d_obs = Taper(obs_rz[nt_start:nt_end],0.05)
                        adj_rz[nt_start:nt_end] = self.Adj(d_syn,d_obs,ntype)
                    adj_rz = self.filt(self.filt(adj_rz,self.par.freq,order=1)[::-1],self.par.freq,order=1)[::-1]
                    adj_rz.astype('f4').tofile(stat[i,0].strip()+'-ADJZ.bin')
                elif nwin[2] == 0:
                    adj_rz.astype('f4').tofile(stat[i,0].strip()+'-ADJZ.bin')
                elif nwin[2] < 0:
                    adj_rz = self.Adj(rz,obs_rz,ntype)
                    adj_rz.astype('f4').tofile(stat[i,0].strip()+'-ADJZ.bin')


    def Error(self,ntype,sufix,iter):
        # Misfit 
        Madj_rx = 0.0
        Madj_ry = 0.0
        Madj_rz = 0.0

        f = open(self.par.FoldIter+'/F-Component'+'.txt', 'w')
        f.close()
        
        
        for k in range(0,self.par.EventNum):
            os.chdir(self.par.InvPar["-WorkDir"])
            os.chdir(self.par.EventName[k])
            # Call Windows Class
            win = windows(self.par)
            
            os.chdir(self.par.EventName[k])

            #Station
            stat = np.genfromtxt('receptors.par',dtype='str',comments='#',delimiter=',')
            n = stat.shape # n[0] number of stations that record the event

            os.chdir('DATA/')
            
            for i in range(0,n[0]):
                rx = np.fromfile(stat[i,0].strip()+'-VX.bin',dtype='f4')
                ry = np.fromfile(stat[i,0].strip()+'-VY.bin',dtype='f4')
                rz = np.fromfile(stat[i,0].strip()+'-VZ.bin',dtype='f4')
                obs_rx = np.fromfile(stat[i,0].strip()+'-VX-OBS.bin',dtype='f4')
                obs_ry = np.fromfile(stat[i,0].strip()+'-VY-OBS.bin',dtype='f4')
                obs_rz = np.fromfile(stat[i,0].strip()+'-VZ-OBS.bin',dtype='f4')



                # Filter Seismograms

                rx = self.filt(Taper(rx,0.05),self.par.freq)
                ry = self.filt(Taper(ry,0.05),self.par.freq)
                rz = self.filt(Taper(rz,0.05),self.par.freq)

                obs_rx = self.filt(Taper(obs_rx,0.05),self.par.freq)
                obs_ry = self.filt(Taper(obs_ry,0.05),self.par.freq)
                obs_rz = self.filt(Taper(obs_rz,0.05),self.par.freq)

                stat_name = stat[i,0]
                # Number of windows [east,north,z]
                nwin = win.stat_win[stat_name]

                # East Component 
                if nwin[0] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    inf_w = win.stat_we[stat_name]
                    for i_nw in range(0,nwin[0]):
                        nt_start = int(np.floor(inf_w[i_nw,0] / self.par.dt))
                        nt_end = int(np.floor(inf_w[i_nw,1] / self.par.dt))
                        d_syn = rx[nt_start:nt_end]
                        d_obs = obs_rx[nt_start:nt_end]
                        Madj_rx += self.Misfit(d_syn,d_obs,ntype)
                elif nwin[0] == 0:
                    Madj_rx += 0.0
                elif nwin[0] < 0:
                    Madj_rx += self.Misfit(rx,obs_rx,ntype)



                # North Component 
                if nwin[1] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    inf_w = win.stat_wn[stat_name]
                    for i_nw in range(0,nwin[1]):
                        nt_start = int(np.floor(inf_w[i_nw,0] / self.par.dt))
                        nt_end = int(np.floor(inf_w[i_nw,1] / self.par.dt))
                        d_syn = ry[nt_start:nt_end]
                        d_obs = obs_ry[nt_start:nt_end]
                        Madj_ry += self.Misfit(d_syn,d_obs,ntype)
                elif nwin[1] == 0:
                    Madj_ry += 0.0
                elif nwin[1] < 0:
                    Madj_ry += self.Misfit(ry,obs_ry,ntype)


                # Z Component 
                if nwin[2] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    inf_w = win.stat_wz[stat_name]
                    for i_nw in range(0,nwin[2]):
                        nt_start = int(np.floor(inf_w[i_nw,0] / self.par.dt))
                        nt_end = int(np.floor(inf_w[i_nw,1] / self.par.dt))
                        d_syn = rz[nt_start:nt_end]
                        d_obs = obs_rz[nt_start:nt_end]
                        Madj_rz += self.Misfit(d_syn,d_obs,ntype)
                elif nwin[2] == 0:
                    Madj_rz += 0.0
                elif nwin[2] < 0:
                    Madj_rz += self.Misfit(rz,obs_rz,ntype)


        f = open(self.par.FoldIter+'/F-Component'+'.txt', 'w')
        f.write('# Misfit ' + ': Ex,Ey,Ez\n')
        f.write('%e,%e,%e\n' % ((Madj_rx),(Madj_ry),(Madj_rz)))
        f.close()


        # Total Misfit
        Mt = (Madj_rx + Madj_ry + Madj_rz) * 0.5  
  
        # if iter > 0:
        #     Mt = (Madj_rx + Madj_ry + Madj_rz) * 0.5
        #     mtest = np.fromfile(self.par.FoldIter +"/m-test.bin",dtype='f4')
        #     mnew = np.fromfile(self.par.FoldIter +"/m-new.bin",dtype='f4')
        #     fnew = np.loadtxt(self.par.FoldIter+'/f-new.txt')
        #     alpha = fnew * 0.1 / (self.par.nodz * self.par.nody * self.par.nodx)
        #     print "ALPHA: ",alpha
        #     var = 1.0 / (0.15* mnew)**2
        #     e_m = alpha * 0.5 * np.linalg.norm((np.sqrt(var) * (mtest - mnew)))**2 
        #     Mt += e_m
        
        f = open(self.par.FoldIter+'/f-'+sufix+'.txt', 'w')
        f.write('%e\n' % Mt)
        f.close()

                    
    def Adj(self,syn,obs,ntype):
        """ ntype options:
        ntype = 1 -> L2 Norm 
        """

        # L2 Norm
        if ntype == 1:
            diff = Taper(syn - obs,0.05)
            adj = diff
        # L2 Norm Frequency
        if ntype == 2:
            adj = FreqExtract(syn,self.par.nt,self.par.dt,self.par.freq) - \
                  FreqExtract(obs,self.par.nt,self.par.dt,self.par.freq)
        if ntype == 3:
            adj = -1.0 * TravelTimeAdj(syn, obs, syn.size, self.par.dt)
        return adj


    def Misfit(self,syn,obs,ntype):
        """ ntype options:
        ntype = 1 -> L2 Norm 
        """
        # L2 Norm
        if ntype == 1:
            wadj = syn - obs
        # L2 Norm Frequency
        if ntype == 2:
            wadj = FreqExtract(syn,self.par.nt,self.par.dt,self.par.freq) - \
                  FreqExtract(obs,self.par.nt,self.par.dt,self.par.freq)
        if ntype == 3:
            wadj = TravelTimeMisfit(syn, obs, syn.size, self.par.dt)
            
        return (np.sum(wadj*wadj))
        

    def filt(self,data,freq,order=4):
        data_filt = self.filter(data,freq[0],freq[-1],1.0/self.par.dt,order)
        return data_filt



def Taper(data,alpha=0.25):
    win = signal.tukey(data.size,alpha)
    return win * data


def TravelTimeAdj(syn, obs, nt, dt):
    # cross correlation traveltime
    # (Tromp et al 2005, eq 45)
    wadj = syn #np.zeros(nt)
    #wadj[1:-1] = (syn[2:] - syn[0:-2])/(2.*dt)
    wadj *= 1./(sum(wadj*wadj)*dt)
    wadj *= -TravelTimeMisfit(syn,obs,nt,dt)
    return wadj

def TravelTimeMisfit(syn, obs, nt, dt):
    cc = abs(np.convolve(obs, np.flipud(syn)))
    return (np.argmax(cc)-nt+1)*dt


def FreqExtract(data,nfd,dt,freq):
        tfd = np.arange(0,nfd)*dt+dt/2.0
        exc_data = data * 0.0
        for ifreq in freq:
                #print ifreq
                fun1 = np.cos(2.0 * np.pi * ifreq * tfd)
                fun2 = np.sin(2.0 * np.pi * ifreq * tfd)
                w0 = np.trapz(data * fun1) * dt
                w1 = np.trapz(data * fun2) * dt
                exc_data += fun1 * w0 + fun2 * w1
        return Taper(exc_data,0.05)    
