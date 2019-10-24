import numpy as np 
import matplotlib.pyplot as plt
import math
from scipy.signal import butter, lfilter
from scipy.signal import freqz,freqz_zpk,sosfilt,zpk2sos
from scipy import signal
import subprocess
import os
import time,sys


## Source time function ##                       
def source(f0,t0,t,dt,t_half):
    alpha = math.pi * math.pi * f0 * f0;
    src = np.exp(-alpha * (t-t0)**2.0)# * np.sqrt(alpha / np.pi);
    #dsrc = 4.0 * alpha *(t - t0) * np.exp(-2.0 * alpha * (t  - t0)**2 );
    dsrc =  np.sqrt(alpha / np.pi) * -2.0 * alpha * (t - t0) * np.exp(-alpha * (t0 - t)**2);
    ddsrc = (1.0 - 2.0 * alpha * (t - t0)**2.0) * np.exp(-alpha * (t - t0)**2.0)

    dsrc_w = np.fft.fftshift(np.fft.fft(dsrc));
    ddsrc_w = np.fft.fftshift(np.fft.fft(ddsrc));

    f0_w = np.fft.fftshift(np.fft.fftfreq(t.shape[0],dt));

    trg = np.array([0.0,1.0,0.0]);
    t_trg = np.array([t0-t_half,t0,t0+t_half])

    plt.close('all')
    plt.figure()
    plt.plot(t,src,'r.-')
    plt.plot(t_trg,trg,'b.-')
    plt.plot(t,dsrc,'b.-')
    plt.title('Gaussian Function',fontsize=9)
    plt.xlabel(r'$t[s]$')
    plt.ylabel(r'$f[t]$')
    plt.grid(True)
    plt.xlim([0,10])
    
    return (dsrc);


def spectrum(src,dt):
    dsrc_w = np.fft.rfft(src)
    f0_w = np.fft.rfftfreq(src.shape[0],dt)
    return [f0_w , dsrc_w]

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    z, p,k = butter(order, [low, high], btype='bandpass',output='zpk')
    return z, p,k

def filter(data, lowcut, highcut, fs, order=5):
    z, p, k = butter_bandpass(lowcut, highcut, fs, order=order)
    sos = zpk2sos(z,p,k)
    y = sosfilt(sos, data)
    return y

def Taper(data,alpha=0.25):
    win = signal.tukey(data.size,alpha)
    return win * data


shot = 'shot5'
command = "ls -l "+ shot +"/DATA/*VX-OBS.bin | awk '{print $9}' | awk -F- '{print $1}' | awk -F/ '{print $3}'"
#Station
stat =subprocess.check_output(command,shell=True).split('\n')

freq = [0.1,2.0]
### Finite Diference
dt = 0.01;
nfd = 600;
t0_fd = 0.0
tfd = np.arange(0,nfd)*dt
fd_fold = shot+'/DATA/'


# Observed Data 
nobs = nfd
dt_obs = dt;
t_obs = np.arange(0,nobs)*dt_obs;
obs_fold = fd_fold



for i in range(0,len(stat)-1,5):

    g_x = Taper(np.fromfile(fd_fold + stat[i] +'-VX.bin',dtype='f4'),0.2)
    g_y = Taper(np.fromfile(fd_fold + stat[i] +'-VY.bin',dtype='f4'),0.2)
    g_z = Taper(np.fromfile(fd_fold + stat[i] +'-VZ.bin',dtype='f4'),0.2)


    fd_y = g_y#np.convolve(g_y,st) * dt
    fd_x = g_x#np.convolve(g_x,st) * dt
    fd_z = g_z#np.convolve(g_z,st) * dt


    # [fw_dft , spec_x] = spectrum(fd_x,dt)
    # [fw_dft , spec_y] = spectrum(fd_y,dt)
    # [fw_dft , spec_z] = spectrum(fd_z,dt)
    
    # fd_y = np.cumsum(np.convolve(g_y,st) * dt) * dt
    # fd_x = np.cumsum(np.convolve(g_x,st) * dt) * dt
    # fd_z = np.cumsum(np.convolve(g_z,st) * dt) * dt


    # fd_y = filter(fd_y, 0.01, 0.05, 1.0/dt, order=4)[n0-1:n0 + nfd-1]
    # fd_x = filter(fd_x, 0.01, 0.05, 1.0/dt, order=4)[n0-1:n0 + nfd-1]
    # fd_z = filter(fd_z, 0.01, 0.05, 1.0/dt, order=4)[n0-1:n0 + nfd-1]

    fd_x = filter(fd_x, freq[-0],freq[-1], 1.0/dt, order=4)
    fd_y = filter(fd_y,freq[-0],freq[-1], 1.0/dt, order=4)
    fd_z = filter(fd_z,freq[-0],freq[-1], 1.0/dt, order=4)

    [fw_dft , spec_x] = spectrum(fd_x,dt)
    [fw_dft , spec_y] = spectrum(fd_y,dt)
    [fw_dft , spec_z] = spectrum(fd_z,dt)

    r_x =  np.fromfile(obs_fold + stat[i] +'-VX-OBS.bin',dtype='f4')
    r_y =  np.fromfile(obs_fold + stat[i] +'-VY-OBS.bin',dtype='f4')
    r_z =  np.fromfile(obs_fold + stat[i] +'-VZ-OBS.bin',dtype='f4')

    r_x =  Taper(r_x,0.2)
    r_y =  Taper(r_y,0.2)
    r_z =  Taper(r_z,0.2)

    # r_y = np.cumsum(r_y) * dt
    # r_x = np.cumsum(r_x) * dt
    # r_z = np.cumsum(r_z) * dt

    # [fw , spec_rx] = spectrum(r_x,dt)
    # [fw , spec_ry] = spectrum(r_y,dt)
    # [fw , spec_rz] = spectrum(r_z,dt)


    r_x =  filter(r_x,freq[-0],freq[-1], 1.0/dt_obs, order=4) 
    r_y =  filter(r_y,freq[-0],freq[-1], 1.0/dt_obs, order=4) 
    r_z =  filter(r_z,freq[-0],freq[-1], 1.0/dt_obs, order=4) 

    [fw , spec_rx] = spectrum(r_x,dt)
    [fw , spec_ry] = spectrum(r_y,dt)
    [fw , spec_rz] = spectrum(r_z,dt)

    # adj_x = r_x[:nfd] - fd_x
    # adj_y = r_y[:nfd] - fd_y
    # adj_z = r_z[:nfd] - fd_z

    # e_x += np.sum(adj_x * adj_x)
    # e_y += np.sum(adj_y * adj_y)
    # e_z += np.sum(adj_z * adj_z) 

    print stat[i]


    plt.figure()
    plt.subplot(9,1,1)
    plt.title(stat[i])
    plt.plot(tfd,fd_x,'r-',label = 'fd')
    plt.plot(t_obs,r_x,'k-',label = 'observed')
    plt.subplot(9,1,2)
    plt.plot(fw_dft,np.abs(spec_x),'r-*',label = 'fd')
    plt.plot(fw,np.abs(spec_rx),'k-*',label = 'observed')
    plt.xlim([0,3])
    plt.subplot(9,1,3)
    plt.plot(fw_dft,(np.angle(spec_x)),'r-*',label = 'fd')
    plt.plot(fw,(np.angle(spec_rx)),'k-*',label = 'observed')
    plt.xlim([0,3])
    #plt.yscale('log')
    plt.legend()
    plt.subplot(9,1,4)    
    plt.plot(tfd,fd_y,'r-')
    plt.plot(t_obs,r_y,'k-')
    plt.subplot(9,1,5)
    plt.plot(fw_dft,np.abs(spec_y),'r-*',label = 'fd')
    plt.plot(fw,np.abs(spec_ry),'k-*',label = 'observed')
    plt.xlim([0,3])
    plt.subplot(9,1,6)
    plt.plot(fw_dft,(np.angle(spec_y)),'r-*',label = 'fd')
    plt.plot(fw,(np.angle(spec_ry)),'k-*',label = 'observed')
    plt.xlim([0,3])
    #plt.yscale('log')
    #plt.xlim([0,50])
    plt.subplot(9,1,7) 
    plt.plot(tfd,fd_z,'r-')   
    plt.plot(t_obs,r_z,'k-')
    plt.subplot(9,1,8)
    plt.plot(fw_dft,np.abs(spec_z),'r-*',label = 'fd')
    plt.plot(fw,np.abs(spec_rz),'k-*',label = 'observed')
    plt.xlim([0,3])
    plt.subplot(9,1,9)
    plt.plot(fw_dft,(np.angle(spec_z)),'r-*',label = 'fd')
    plt.plot(fw,(np.angle(spec_rz)),'k-*',label = 'observed')
    plt.xlim([0,3])
    #plt.yscale('log')
    #plt.xlim([0,50])


plt.show()

# print "e_x,e_y,e_z",e_x,e_y,e_z

 
