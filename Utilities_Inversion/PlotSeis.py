import numpy as np 
import matplotlib.pyplot as plt
import math
import sys
from scipy import signal


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
        return exc_data


sys.path.append("inversion/")
from Tools.UtilFunc import butter_bandpass_filter,PlotNorm,butter_lowpass_filter

#Parameters Modeling 
info = np.genfromtxt('ParForward.par',dtype='str',comments='//',delimiter=' ')

### Finite Diference
dt = float(info[10,1])
nfd = int(float(info[9,1])/dt)
tfd = np.arange(0,nfd)*dt+dt/2.0
fold = "shot19/"
filter = butter_bandpass_filter
freq = 2.0

#print dt,nfd,tfd


#Station
stat = np.genfromtxt(fold + 'receptors.par',dtype='str',comments='#',delimiter=',')
nstat = stat.shape[0] # n[0] number of stations that record the event

#print stat


#Source
source = np.genfromtxt(fold + 'sources.par',dtype='str',comments='#',delimiter=',')
nsour = stat.shape # n[0] number of stations that record the event

xsrc = float(source[0])
ysrc = float(source[1])
zsrc = float(source[2])

ifreq = [0.1, 0.9,3.0]

print ifreq


dist = np.zeros((nstat))

for i in range(0,nstat):
        xstat = float(stat[i,1])
        ystat = float(stat[i,2])
        zstat = float(stat[i,3])

        dist[i] = np.sqrt( (xstat - xsrc)**2 +\
                        (ystat - ysrc)**2 +\
                        (zstat - zsrc)**2)

order = np.argsort(dist)
#print order
#print dist[order]         
        

#print source

for i in range(0,nstat):
	rx = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VX-M0.bin',dtype='f4')
	ry = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VY-M0.bin',dtype='f4')
	rz = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VZ-M0.bin',dtype='f4')


        rx_i = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VX.bin',dtype='f4')
	ry_i = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VY.bin',dtype='f4')
	rz_i = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VZ.bin',dtype='f4')
        

        #rx = FreqExtract(rx,nfd,dt,ifreq)
        #ry = FreqExtract(ry,nfd,dt,ifreq)
        #rz = FreqExtract(rz,nfd,dt,ifreq)

        obs_rx = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VX-OBS.bin',dtype='f4')
	obs_ry = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VY-OBS.bin',dtype='f4')
	obs_rz = np.fromfile(fold +'DATA/' + stat[i,0].strip()+'-VZ-OBS.bin',dtype='f4')


        # adj_x1 = rx - obs_rx

        # adj_x1_f = filter(filter(adj_x1,ifreq[-1],1/dt,5)[::-1],ifreq[-1],1/dt,5)[::-1]

        rx = filter(rx,ifreq[0],ifreq[-1],1/dt,4)
        ry = filter(ry,ifreq[0],ifreq[-1],1/dt,4)
        rz = filter(rz,ifreq[0],ifreq[-1],1/dt,4)

        rx_i = filter(rx_i,ifreq[0],ifreq[-1],1/dt,4)
        ry_i = filter(ry_i,ifreq[0],ifreq[-1],1/dt,4)
        rz_i = filter(rz_i,ifreq[0],ifreq[-1],1/dt,4)


        obs_rx = filter(obs_rx,ifreq[0],ifreq[-1],1/dt,4)
        obs_ry = filter(obs_ry,ifreq[0],ifreq[-1],1/dt,4)
        obs_rz = filter(obs_rz,ifreq[0],ifreq[-1],1/dt,4)

        adj_x = dt * (np.cumsum(rx) - np.cumsum(obs_rx))
        adj_y = dt * (np.cumsum(ry) - np.cumsum(obs_ry))
        adj_z = dt * (np.cumsum(rz) - np.cumsum(obs_rz))

        #obs_rx = FreqExtract(obs_rx,nfd,dt,ifreq)
        #obs_ry = FreqExtract(obs_ry,nfd,dt,ifreq)
        #obs_rz = FreqExtract(obs_rz,nfd,dt,ifreq)


        #rx *= 1.0 / np.max(np.abs(obs_rx))
        #ry *= 1.0 / np.max(np.abs(obs_ry))
        #rz *= 1.0 / np.max(np.abs(obs_rz))

        #obs_rx *= 1.0 / np.max(np.abs(obs_rx))
        #obs_ry *= 1.0 / np.max(np.abs(obs_ry))
        #obs_rz *= 1.0 / np.max(np.abs(obs_rz))

        # adj_rx = np.fromfile(fold +'DATA/' +  stat[i,0].strip()+'-ADJX.bin',dtype='f4')
	# adj_ry = np.fromfile(fold +'DATA/' +  stat[i,0].strip()+'-ADJY.bin',dtype='f4')
	# adj_rz = np.fromfile(fold +'DATA/' +  stat[i,0].strip()+'-ADJZ.bin',dtype='f4')

        #adj_rx *= 1.0 / np.max(np.abs(adj_rx))
        #adj_ry *= 1.0 / np.max(np.abs(adj_ry))
        #adj_rz *= 1.0 / np.max(np.abs(adj_rz))


        plt.figure(1)
        plt.subplot(3,1,1)
	plt.plot(tfd,rx,'b-')
        plt.plot(tfd,rx_i,'r-')
	plt.plot(tfd,obs_rx,'k-') 
	plt.text(0.0,0.0,stat[i,0].strip()+'-E',fontsize=12)

        plt.subplot(3,1,2)       
	plt.plot(tfd,ry,'b-')
        plt.plot(tfd,ry_i,'r-') 
	plt.plot(tfd,obs_ry,'k-')
	plt.text(0.0,0.0,stat[i,0].strip()+'-N',fontsize=12)

        plt.subplot(3,1,3)         
	plt.plot(tfd,rz,'b-')
        plt.plot(tfd,rz_i,'r-')
	plt.plot(tfd,obs_rz,'k-')
	plt.text(0.0,0.0,stat[i,0].strip()+'-Z',fontsize=12)



        #adj_f = FreqExtract(rx - obs_rx,nfd,dt,ifreq)
        #adj = filter(adj_f[::-1],ifreq[-1],1/dt,5)[::-1]
        
        # plt.figure(2)
        # plt.subplot(3,1,1)
        # #plt.plot(tfd,rx,'b-')
	# #plt.plot(tfd,adj_x,'k-')
        # plt.plot(tfd,adj_rx,'ro')
	# plt.text(0.0,0.0,stat[i,0].strip()+'-E',fontsize=12)

        # plt.subplot(3,1,2) 
        # #plt.plot(tfd,adj_y,'k-')
	# plt.plot(tfd,adj_ry,'r-')
	# plt.text(0.0,0.0,stat[i,0].strip()+'-N',fontsize=12)

        # plt.subplot(3,1,3) 
        # #plt.plot(tfd,adj_z,'k-')
	# plt.plot(tfd,adj_rz,'r-')
	# plt.text(0.0,0.0,stat[i,0].strip()+'-Z',fontsize=12)
        


        # plt.figure(1)            
	# plt.plot(tfd,rx + 2.5 * order[i],'r-')
	# plt.plot(tfd,obs_rx + 2.5 * order[i],'k-') 
	# plt.text(0.1,2.5 * order[i] + .1,stat[i,0].strip()+'-E',fontsize=8)

        # plt.figure(2)            
	# plt.plot(tfd,ry + 2.5 * order[i],'r-')
	# plt.plot(tfd,obs_ry + 2.5 * order[i],'k-')
	# plt.text(0.1,2.5 * order[i] + .1,stat[i,0].strip()+'-N',fontsize=8)

        # plt.figure(3)            
	# plt.plot(tfd,rz + 2.5 * order[i],'r-')
	# plt.plot(tfd,obs_rz + 2.5 *  order[i],'k-')
	# plt.text(0.1,2.5 * order[i] +.1,stat[i,0].strip()+'-Z',fontsize=8)

        # plt.figure(4)            
	# plt.plot(tfd,adj_rx + 2.5 * order[i],'b-')
        # plt.plot(tfd,rx - obs_rx + 2.5 * order[i],'r-')
	# plt.text(0.1,2.5 * order[i] + .1,stat[i,0].strip()+'-ADJE',fontsize=8)

        # plt.figure(5)            
	# plt.plot(tfd,adj_ry + 2.5 * order[i],'b-')
        # plt.plot(tfd,ry - obs_ry + 2.5 * order[i],'r-')
	# plt.text(0.1,2.5 * order[i] + .1,stat[i,0].strip()+'-ADJN',fontsize=8)

        # plt.figure(6)            
	# plt.plot(tfd,adj_rz + 2.5 * order[i],'b-')
        # plt.plot(tfd,rz - obs_rz + 2.5 * order[i],'r-')
	# plt.text(0.1,2.5 * order[i] +.1,stat[i,0].strip()+'-ADJZ',fontsize=8)
        plt.show()
        plt.close('all')




        
	# plt.subplot(3,1,2)    
	# plt.plot(tfd,ry,'r-',label = stat[i,0].strip()+'-N')
        # plt.plot(tfd,obs_ry,'k-')
 	# plt.legend(fontsize=8,loc=2)
	# plt.subplot(3,1,3) 
	# plt.plot(tfd,rz,'r-',label = stat[i,0].strip()+'-Z')
        # plt.plot(tfd,obs_rz,'k-')
 	# plt.legend(fontsize=8,loc=2)
        # plt.show()




