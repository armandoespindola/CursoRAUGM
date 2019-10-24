import numpy as np 
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from scipy.signal import freqz,freqz_zpk,sosfilt,zpk2sos
import obspy
#from obspy.signal.tf_misfit import plot_tf_misfits


 ### Show Images Kernel
def Simshow(K,cx,cy,cz,sx,sy,sz,title):
    plt.figure()
    Dy = cy[1]-cy[0]
    Dx = cx[1]-cx[0]
    Dz = cz[0]-cz[1]
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=3)
    plt.subplot(3,1,1)
    plt.title(title+'-ZY',fontsize=5)
    plt.set_cmap('seismic')
    #plt.contourf(K[:,:,sx],levels=lev0)
    plt.imshow(K[:,:,sx],origin='lower',aspect='equal')
    plt.xticks(cy[::30]/Dy,[str(cy[i]/1000.0) for i in range(0,cy.size,30)])
    plt.yticks(cz[::20]/Dz,[str((cz[0]-cz[i])/1000.0) for i in range(0,cz.size,20)])
    plt.colorbar(format = '%.2E')
    

    plt.subplot(3,1,2)
    plt.title(title+'-ZX',fontsize=5)
    plt.set_cmap('seismic')
    #plt.contourf(K[:,sy,:],levels=lev1)
    plt.imshow(K[:,sy,:],origin='lower')
    plt.xticks(cx[::30]/Dx,[str(cx[i]/1000.0) for i in range(0,cx.size,30)])
    plt.yticks(cz[::20]/Dz,[str((cz[0]-cz[i])/1000.0) for i in range(0,cz.size,20)])
    plt.colorbar(format = '%.2E')
 

    plt.subplot(3,1,3)
    plt.title(title+'-YX',fontsize=5)
    plt.set_cmap('seismic')
    #plt.contourf(K[sz,:,:],levels=lev2)
    plt.imshow(K[sz,:,:],origin='lower')
    plt.xticks(cy[::30]/Dy,[str(cy[i]/1000.0) for i in range(0,cy.size,30)])
    plt.yticks(cx[::30]/Dx,[str(cx[i]/1000.0) for i in range(0,cx.size,30)])
    plt.colorbar(format = '%.2E')
    plt.subplots_adjust(hspace = 0.4)


### Source time function ###                       

def source(f0,t0,time,dt,t_half,Type):
    nt = round(time/dt)
    t = np.arange(0,nt) * dt + dt/2.0;
    alpha = math.pi * math.pi * f0 * f0;
    src = np.exp(-alpha * (t-t0)**2.0);            
    dsrc = 4.0 * alpha *(t - t0) * np.exp(-2.0 * alpha * (t  - t0)**2 );
    #dsrc = (t - t0) * np.exp(-alpha * (t0 - t)**2);
    ddsrc = (1.0 - 2.0 * alpha * (t - t0)**2.0) * np.exp(-alpha * (t - t0)**2.0)

    if Type==1:
        src = src
    elif Type==2:
        src = dsrc
    elif Type==3:
        src = ddsrc
    else:
        print "this is not an option Source Type"
    
    trg = np.array([0.0,1.0,0.0]);
    t_trg = np.array([t0-t_half,t0,t0+t_half])

    plt.close('all')
    plt.figure()
    plt.plot(t,src,'r.-')
    plt.title('Source Time',fontsize=9)
    plt.xlabel(r'$t[s]$')
    plt.grid(True)
    plt.savefig("SourceTime.pdf",dpi=300,bbox_inches='tight')
    
    return (src);


### Butterworth Bandpass ###

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    z, p,k = butter(order, [low, high], btype='bandpass',output='zpk')
    return z, p,k

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    z, p, k = butter_bandpass(lowcut, highcut, fs, order=order)
    sos = zpk2sos(z,p,k)
    y = sosfilt(sos, data)
    return y

### Butterworth LowPass ###

def butter_lowpass(lowcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    z, p,k = butter(order, [low], btype='lowpass',output='zpk')
    return z, p,k

def butter_lowpass_filter(data, lowcut, fs, order=5):
    z, p, k = butter_lowpass(lowcut, fs, order=order)
    sos = zpk2sos(z,p,k)
    y = sosfilt(sos, data)
    return y


def PlotNorm(rx,obs_rx,ry,obs_ry,rz,obs_rz,adj_rx,adj_ry,adj_rz,title,par):
    dt = par.dt;
    nfd = par.Tsim / dt;
    tfd = np.arange(0,nfd)*dt+dt/2.0;

    plt.figure()
    plt.subplot(3,1,1)
    plt.plot(tfd,rx,'b-',label = 'Syn-E')
    plt.plot(tfd,obs_rx,'r-',label = 'Obs-E')
    plt.legend(fontsize=8,loc=2)
    plt.subplot(3,1,2)
    plt.plot(tfd,ry,'b-',label = 'Syn-N')
    plt.plot(tfd,obs_ry,'r-',label = 'Obs-N')
    plt.legend(fontsize=8,loc=2)
    plt.subplot(3,1,3) 
    plt.plot(tfd,rz,'b-',label ='Syn-Z')
    plt.plot(tfd,obs_rz,'r-',label ='Obs-Z')
    plt.legend(fontsize=8,loc=2)
    plt.savefig(title+".pdf",dpi=300,bbox_inches='tight');
     
    plt.figure()
    plt.subplot(3,1,1)
    plt.plot(tfd,adj_rx,'k-',label = 'Adj-E')
    plt.legend(fontsize=8,loc=2)
    plt.subplot(3,1,2)    
    plt.plot(tfd,adj_ry,'k-',label = 'Adj-N')
    plt.legend(fontsize=8,loc=2)
    plt.subplot(3,1,3) 
    plt.plot(tfd,adj_rz,'k-',label ='Adj-Z')
    plt.legend(fontsize=8,loc=2)
    plt.savefig(title+"-ADJ.pdf",dpi=300,bbox_inches='tight');

    plt.close('all')


#def FreqMisfit(rx,obs_rx,ry,obs_ry,rz,obs_rz,freqmin,freqmax,title,par):

#    dt = par.dt;

#    plt.figure()
#    plot_tf_misfits(rx, obs_rx, dt=dt, fmin=freqmin, fmax=freqmax, show=False)
#    plt.savefig(title+"FreqM-x.pdf",dpi=300,bbox_inches='tight');

#    plt.figure()
#    plot_tf_misfits(ry, obs_ry, dt=dt, fmin=freqmin, fmax=freqmax, show=False)
#    plt.savefig(title+"FreqM-y.pdf",dpi=300,bbox_inches='tight');

#    plt.figure()
#    plot_tf_misfits(rz, obs_rz, dt=dt, fmin=freqmin, fmax=freqmax, show=False)
#    plt.savefig(title+"FreqM-z.pdf",dpi=300,bbox_inches='tight');

#    plt.close('all')



def polyfitB(x,f):
    # parabolic fit
    i = np.argmin(f)
    p = np.polyfit(x[i-1:i+2], f[i-1:i+2], 2)

    if p[0] > 0:
        return -p[1]/(2*p[0])
    else:
        print -1
        raise Exception()

def angle(x,y):
    xy = dot(x,y)
    xx = dot(x,x)
    yy = dot(y,y)
    b = xy/((xx*yy)**0.5)
    if b >= 0.99999:
        b = 0.99990
    #print "xx,yy,xy,arg: ",xx,yy,xy,b
    return np.arccos(b)


def dot(x,y):
    return np.dot(
        np.squeeze(x),
        np.squeeze(y))



def backtrack2(f0, g0, x1, f1, b1=0.1, b2=0.5):
    """ Safeguarded parabolic backtrack
        Ryan Modrak 2016;
    """
    # parabolic backtrack
    x2 = -g0*x1**2/(2*(f1-f0-g0*x1))

    

    # apply safeguards
    if x2 > b2*x1:
        x2 = b2*x1
    elif x2 < b1*x1:
        x2 = b1*x1
    return x2


def count_zeros(a):
    """ Counts number of zeros in a list or array
    """
    return sum(np.array(a)==0)

def save(name,var):
    np.savetxt(name, [var], fmt='%10.4e')

def fft(data,dt):
    fdata = np.fft.fftshift(np.fft.fft(data))
    f0_w = np.fft.fftshift(np.fft.fftfreq(data.shape[0],dt))
    return [fdata , f0_w]


def bin2sac(data,dt,NameEvent,Component,StatName):
    # Create Obspy Trace
    tr = obspy.core.trace.Trace(np.float64(data))
    tr.stats.network = NameEvent
    tr.stats.delta = dt
    tr.stats.samplin_rate = 1.0 / dt
    tr.stats.station = StatName+ "-" + Component
    tr.stats.channel = Component
    tr.stats.sampling_rate = 1 / dt
    stream = obspy.core.stream.Stream(traces=[tr])
    return stream

def CloneBin2Sac(obs,clone):
    syn_data = obs.copy()
    syn_data[0].data = clone
    return syn_data
    
    
