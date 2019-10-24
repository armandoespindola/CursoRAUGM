import numpy as np
import os
from inversion.Tools.UtilFunc import butter_bandpass_filter,PlotNorm,butter_lowpass_filter,bin2sac,CloneBin2Sac
import math
from inversion.Windows.windows import windows
from scipy import signal
import time
import obspy

class AdjointSource:
    def __init__(self,par,veldisp,syn_sac,obs_sac):
        self.par = par
        print '## freq: ',self.par.freq
        self.filter = butter_bandpass_filter
        self.fw = np.fft.rfftfreq(self.par.nt,self.par.dt)
        self.veldisp = veldisp
        self.syn_sac = syn_sac
        self.obs_sac = obs_sac
        self.nwindows = 0
        #print self.fw
        
    def ReadDir(self,i):
        os.chdir(self.par.InvPar["-WorkDir"])
        os.chdir(self.par.EventName[i] + "/DATA")
        

    def AdjSrc(self,ntype,iter):
        freq = self.par.freq[-1]
        print freq
        
        for k in range(0,self.par.EventNum):
            os.chdir(self.par.InvPar["-WorkDir"])
            os.chdir(self.par.EventName[k])
            print self.par.EventName[k]

            #Station
            stat = np.genfromtxt('receptors.par',dtype='str',comments='#',delimiter=',')
            n = stat.shape # n[0] number of stations that record the event

            os.chdir('DATA/')
            
            for i in range(0,n[0]):
                #print stat[i,0]

                name_stat = stat[i,0].strip()
                rx_v,ry_v,rz_v,obs_rx_v,obs_ry_v,obs_rz_v = ReadData(name_stat,self.par.dt\
                                                                     ,self.syn_sac,self.obs_sac)

                rx = ProcessData(rx_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                ry = ProcessData(ry_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                rz = ProcessData(rz_v,self.veldisp,0.05,self.par.freq,zerophase=False)

                obs_rx = ProcessData(obs_rx_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                obs_ry = ProcessData(obs_ry_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                obs_rz = ProcessData(obs_rz_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                

                # Initialize Adjoint Variables
                adj_rx = rx[0].data * 0.0
                adj_ry = ry[0].data * 0.0
                adj_rz = rz[0].data * 0.0




                # if iter>=0:
                # # Compute Windows
                #     if self.obs_sac:
                #         self.win.CreateInfoStatEvent(obs_rx,0.0)

              
                #     self.win.CreateWindows(obs_rx,rx,filename=stat[i,0] + 'VX',plot=False)
                #     winx = self.win.windows
                #     self.win.CreateWindows(obs_ry,ry,filename=stat[i,0] + 'VY',plot=False)
                #     winy = self.win.windows
                #     self.win.CreateWindows(obs_rz,rz,filename=stat[i,0] + 'VZ',plot=False)
                #     winz = self.win.windows

                # else:
                #     self.win.LoadWindows(obs_rx,rx,filename=stat[i,0] + 'VX')
                #     winx = self.win.windows
                #     self.win.LoadWindows(obs_ry,ry,filename=stat[i,0] + 'VY')
                #     winy = self.win.windows
                #     self.win.LoadWindows(obs_rz,rz,filename=stat[i,0] + 'VZ')
                #     winz = self.win.windows
                    
                        
                #nwin = [len(winx),len(winy),len(winz)]
                nwin = [-1,-1,-1]
                

                rx = rx[0].data
                ry = ry[0].data
                rz = rz[0].data

                obs_rx = obs_rx[0].data
                obs_ry = obs_ry[0].data
                obs_rz = obs_rz[0].data

                # East Component
                if nwin[0] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    for i_nw in range(0,nwin[0]):
                        nt_start = winx[i_nw].left
                        nt_end = winx[i_nw].right
                        d_syn = Taper(rx[nt_start:nt_end],0.1)
                        d_obs = Taper(obs_rx[nt_start:nt_end],0.1)
                        adj_rx[nt_start:nt_end] = self.Adj(d_syn,d_obs,ntype)
                    adj_rx = self.filt(self.filt(adj_rx,self.par.freq,order=2)[::-1],self.par.freq,order=2)[::-1]
                    adj_rx.astype('f4').tofile(stat[i,0].strip()+'-ADJX.bin')
                    
                elif nwin[0] < 0:
                    adj_rx = self.Adj(rx,obs_rx,ntype)
                    adj_rx.astype('f4').tofile(stat[i,0].strip()+'-ADJX.bin')
                elif nwin[0] == 0:
                    adj_rx.astype('f4').tofile(stat[i,0].strip()+'-ADJX.bin')

                # North Component
                if nwin[1] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    for i_nw in range(0,nwin[1]):
                        nt_start = winy[i_nw].left
                        nt_end = winy[i_nw].right
                        d_syn = Taper(ry[nt_start:nt_end],0.1)
                        d_obs = Taper(obs_ry[nt_start:nt_end],0.1)
                        adj_ry[nt_start:nt_end] = self.Adj(d_syn,d_obs,ntype)
                    adj_ry = self.filt(self.filt(adj_ry,self.par.freq,order=2)[::-1],self.par.freq,order=2)[::-1]
                    adj_ry.astype('f4').tofile(stat[i,0].strip()+'-ADJY.bin')
                elif nwin[1] < 0:
                    adj_ry = self.Adj(ry,obs_ry,ntype)
                    adj_ry.astype('f4').tofile(stat[i,0].strip()+'-ADJY.bin')
                elif nwin[1] == 0:
                    adj_ry[:] = 0.0
                    adj_ry.astype('f4').tofile(stat[i,0].strip()+'-ADJY.bin')


                # Z Component
                if nwin[2] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    for i_nw in range(0,nwin[2]):
                        nt_start = winz[i_nw].left
                        nt_end = winz[i_nw].right
                        d_syn = Taper(rz[nt_start:nt_end],0.1)
                        d_obs = Taper(obs_rz[nt_start:nt_end],0.1)
                        adj_rz[nt_start:nt_end] = self.Adj(d_syn,d_obs,ntype)
                    adj_rz = self.filt(self.filt(adj_rz,self.par.freq,order=2)[::-1],self.par.freq,order=2)[::-1]
                    adj_rz.astype('f4').tofile(stat[i,0].strip()+'-ADJZ.bin')
            
                elif nwin[2] < 0:
                    adj_rz = self.Adj(rz,obs_rz,ntype)
                    adj_rz.astype('f4').tofile(stat[i,0].strip()+'-ADJZ.bin')
                elif nwin[2] == 0:
                    adj_rz[:] = 0.0
                    adj_rz.astype('f4').tofile(stat[i,0].strip()+'-ADJZ.bin')
                    
            

    def Error(self,ntype,sufix,iter):

        Madj_rx = 0.0 
        Madj_ry = 0.0 
        Madj_rz = 0.0
        self.nwindows = 0
        
        for k in range(0,self.par.EventNum):
            os.chdir(self.par.InvPar["-WorkDir"])
            os.chdir(self.par.EventName[k])

            Edx = 0.0 
            Edy = 0.0 
            Edz = 0.0
  
            #Station
            stat = np.genfromtxt('receptors.par',dtype='str',comments='#',delimiter=',')
            n = stat.shape # n[0] number of stations that record the event

            os.chdir('DATA/')

            for i in range(0,n[0]):
                print stat[i,0]

                name_stat = stat[i,0].strip()
                rx_v,ry_v,rz_v,obs_rx_v,obs_ry_v,obs_rz_v = ReadData(name_stat,self.par.dt\
                                                                     ,self.syn_sac,self.obs_sac)

                rx = ProcessData(rx_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                ry = ProcessData(ry_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                rz = ProcessData(rz_v,self.veldisp,0.05,self.par.freq,zerophase=False)

                obs_rx = ProcessData(obs_rx_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                obs_ry = ProcessData(obs_ry_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                obs_rz = ProcessData(obs_rz_v,self.veldisp,0.05,self.par.freq,zerophase=False)
                

                # Initialize Adjoint Variables
                adj_rx = rx[0].data * 0.0
                adj_ry = ry[0].data * 0.0
                adj_rz = rz[0].data * 0.0


                
                # self.win.LoadWindows(obs_rx,rx,filename=stat[i,0] + 'VX')
                # winx = self.win.windows
                # self.win.LoadWindows(obs_ry,ry,filename=stat[i,0] + 'VY')
                # winy = self.win.windows
                # self.win.LoadWindows(obs_rz,rz,filename=stat[i,0] + 'VZ')
                # winz = self.win.windows

                # nwin = [len(winx),len(winy),len(winz)]
                # self.nwindows += len(winx) + len(winy) + len(winz)
                nwin = [-1,-1,-1]


                rx = rx[0].data
                ry = ry[0].data
                rz = rz[0].data

                obs_rx = obs_rx[0].data
                obs_ry = obs_ry[0].data
                obs_rz = obs_rz[0].data

                
                # East Component
                if nwin[0] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    for i_nw in range(0,nwin[0]):
                        nt_start =  winx[i_nw].left
                        nt_end =  winx[i_nw].right
                        d_syn = Taper(rx[nt_start:nt_end],0.1)
                        d_obs = Taper(obs_rx[nt_start:nt_end],0.1)
                        Edx += self.Misfit(d_syn,d_obs,ntype) / nwin[0]
                elif nwin[0] < 0:
                    Edx += self.Misfit(rx,obs_rx,ntype)
                elif nwin[0] == 0:
                    Edx += 0.0
                    
                    
                # North Component
                if nwin[1] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    for i_nw in range(0,nwin[1]):
                        nt_start =  winy[i_nw].left
                        nt_end =  winy[i_nw].right
                        d_syn = Taper(ry[nt_start:nt_end],0.1)
                        d_obs = Taper(obs_ry[nt_start:nt_end],0.1)
                        Edy += self.Misfit(d_syn,d_obs,ntype) / nwin[1]
                elif nwin[1] < 0:
                    Edy += self.Misfit(ry,obs_ry,ntype)
                elif nwin[1] == 0:
                    Edy += 0.0
                    

                # Z Component
                if nwin[2] > 0 :
                    # Infomation Windows [win_start,win_end,Tshift,CC,dlnA]
                    for i_nw in range(0,nwin[2]):
                        nt_start =  winz[i_nw].left
                        nt_end =  winz[i_nw].right
                        d_syn = Taper(rz[nt_start:nt_end],0.1)
                        d_obs = Taper(obs_rz[nt_start:nt_end],0.1)
                        Edz += self.Misfit(d_syn,d_obs,ntype) / nwin[2]
                elif nwin[2] < 0:
                    Edz += self.Misfit(rz,obs_rz,ntype)
                elif nwin[2] == 0:
                    Edz += 0.0
                    
                    


            #np.savetxt('../FEvent.txt',(Edx + Edy + Edz))

            Madj_rx += Edx
            Madj_ry += Edy 
            Madj_rz += Edz
            

        # Total Misfit
        M = (Madj_rx + Madj_ry + Madj_rz)

        #np.savetxt(self.par.FoldIter+'/fw.txt',M)

        Mt = M
  
        
        f = open(self.par.FoldIter+'/f-'+sufix+'.txt', 'w')
        f.write('%e\n' % Mt)
        f.close()

        f = open(self.par.FoldIter+'/Nwin.txt', 'w')
        f.write('%d\n' % self.nwindows)
        f.close()

                    
    def Adj(self,syn,obs,ntype):
        """ ntype options:
        ntype = 1 -> L2 Norm 
        """

        # L2 Norm
        if ntype == 1:
           adj = syn - obs
            
        if ntype == 3:
            adj = -1.0 * TravelTimeAdj(syn, obs, syn.size, self.par.dt)
        return adj


    def Misfit(self,syn,obs,ntype):
        """ ntype options:
        ntype = 1 -> L2 Norm 
        """
        # L2 Norm
        if ntype == 1:
            wadj = np.sum((syn - obs)**2.0) * self.par.dt
            
        if ntype == 3:
            wadj = TravelTimeMisfit(syn, obs, syn.size, self.par.dt)
            wadj = 0.5 * wadj * wadj
            
        return (wadj)
        

    def filt(self,data,freq,order=4):
        data_filt = self.filter(data,freq[0],freq[-1],1.0/self.par.dt,order)
        #data = self.filter(data_filt,freq[0],freq[-1],1.0/self.par.dt,order)[::-1]
        #data_filt = data
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

def TravelTimeFilt(syn, obs, nt, dt,freq):
    syn_f = filto(syn,dt,freq)
    obs_f = filto(obs,dt,freq)
    cc = abs(np.convolve(obs_f, np.flipud(syn_f)))
    return (np.argmax(cc)-nt+1)*dt


def FreqExtract(data):
    dataf = np.fft.rfft(data)
    phase = (np.angle(dataf))
    amp = np.abs(dataf)
    return [amp,phase]    


def filto(data,dt,freq,order=4):
    filterr = butter_bandpass_filter
    data_filt = filterr(data,freq[0],freq[-1],1.0/dt)
    #data = filterr(data_filt,freq[0],freq[-1],1.0/dt)[::-1]
    #data_filt = data
    return data_filt


def ReadData(name,dt,syn_sac,obs_sac):
        
    if (obs_sac):
        obs_x = obspy.read(name + '-VX-OBS.sac')
        obs_y = obspy.read(name + '-VY-OBS.sac')
        obs_z = obspy.read(name + '-VZ-OBS.sac')
        
    else:
        obs_x = np.fromfile(name + '-VX-OBS.bin',dtype='f4')
        obs_y = np.fromfile(name + '-VY-OBS.bin',dtype='f4')
        obs_z = np.fromfile(name + '-VZ-OBS.bin',dtype='f4')
        obs_x = bin2sac(obs_x,dt,"Inversion","VX",name)
        obs_y = bin2sac(obs_y,dt,"Inversion","VY",name)
        obs_z = bin2sac(obs_z,dt,"Inversion","VZ",name)

    if not syn_sac:
        syn_x = np.fromfile(name + '-VX.bin',dtype='f4')
        syn_x = CloneBin2Sac(obs_x,syn_x)
        
        syn_y = np.fromfile(name + '-VY.bin',dtype='f4')
        syn_y = CloneBin2Sac(obs_y,syn_y)
        
        syn_z = np.fromfile(name + '-VZ.bin',dtype='f4')
        syn_z = CloneBin2Sac(obs_z,syn_z)

    
        
    return [syn_x,syn_y,syn_z,obs_x,obs_y,obs_z]

def ProcessData(dat,veldisp,taper,freq,zerophase=False):
    data = dat
    if veldisp == 'disp' :
        data[0].data = np.cumsum(data[0].data) * (1.0 / data[0].stats.sampling_rate)
        data.detrend("linear")
        data.taper(max_percentage=taper, type="hann")
        data.filter("bandpass", freqmin=freq[0], freqmax=freq[-1],
                    corners=4, zerophase=zerophase)
    elif veldisp == 'vel':
        data.detrend("linear")
        data.taper(max_percentage=taper, type="hann")
        data.filter("bandpass", freqmin=freq[0], freqmax=freq[-1],
                    corners=4, zerophase=zerophase)

    return data
        
        

