import numpy as np
import os
import shutil
from scipy import signal
import matplotlib
import matplotlib.pyplot as plt
import subprocess

from inversion.Tools.UtilFunc import Simshow

class ProcKernel:
    def __init__(self,par):
        self.par = par
        self.cx = np.arange(0,self.par.nodx) * self.par.dx
        self.cy = np.arange(0,self.par.nody) * self.par.dy
        self.cz = np.arange(self.par.nodz-1,-1,-1) * self.par.dz
        mnew = np.fromfile(self.par.FoldIter +"/m-new.bin",dtype='f4')
        nstride = self.par.nodx * self.par.nody * self.par.nodz

        self.rho0 = mnew[0:nstride]
        self.vp0 = mnew[nstride:2 * nstride]
        self.vs0 = mnew[2 * nstride:3 * nstride]
        self.max_g = 1.0
        self.minvs = np.min(self.vs0)
        
        
    def SumKernel(self,sufix,iter):
        os.chdir(self.par.InvPar["-WorkDir"]);
        Kr = np.zeros(int(self.par.ModPar["-nx"]) * int(self.par.ModPar["-ny"]) * int(self.par.ModPar["-nz"]))
        Kvp = np.zeros(int(self.par.ModPar["-nx"]) * int(self.par.ModPar["-ny"]) * int(self.par.ModPar["-nz"]))
        Kvs = np.zeros(int(self.par.ModPar["-nx"]) * int(self.par.ModPar["-ny"]) * int(self.par.ModPar["-nz"]))         
        
        for i in range(0,self.par.EventNum):
            os.chdir(self.par.EventName[i] + "/DATA")
            Kr = Kr + np.fromfile("KRHO.bin",dtype='f4')
            #subprocess.call("rm -f KRHO.bin",shell=True)
            Kvp = Kvp + np.fromfile("KVP.bin",dtype='f4')
            #subprocess.call("rm -f KVP.bin",shell=True)
            Kvs = Kvs + np.fromfile("KVS.bin",dtype='f4')
            #subprocess.call("rm -f KVS.bin",shell=True)
            os.chdir(self.par.WorkDir)
        
        os.chdir("IterInversion")
        Kr.astype('f4').tofile("Kr"+sufix+".bin")
        Kvp.astype('f4').tofile("Kvp"+sufix+".bin")
        Kvs.astype('f4').tofile("Kvs"+sufix+".bin")
        os.chdir(self.par.InvPar["-WorkDir"])

        K = np.zeros(int(self.par.ModPar["-nx"]) * int(self.par.ModPar["-ny"]) * int(self.par.ModPar["-nz"]))
        os.chdir(self.par.InvPar["-WorkDir"])
            
        for i in range(0,self.par.EventNum):
            os.chdir(self.par.EventName[i] + "/DATA")
            K = K + np.abs(np.fromfile("PcondA.bin",dtype='f4')) #+ np.fromfile("PcondB.bin",dtype='f4')
            #K = self.TaperSourRec(np.fromfile("PcondB.bin",dtype='f4'))
            #subprocess.call("rm -f PcondA.bin",shell=True)
            #subprocess.call("rm -f PcondB.bin",shell=True)
            os.chdir("../../")

        os.chdir("IterInversion")
        K.astype('f4').tofile("PcondA"+sufix+".bin")
        os.chdir(self.par.InvPar["-WorkDir"])
        
        
    def Ksmoth(self,iter,sufix,mlen1,mlen2):
        # Smoothing Size
        # D1 Kernel --> mlen1
        # D2 Precond --> mlen2

        os.chdir(self.par.InvPar["-WorkDir"])        

        D1 = KernelSmoothing(mlen1,self.par)
        D2 = KernelSmoothing(mlen2,self.par)

        #print minlength, self.par.f0
        
        #Reading Kernels
        os.chdir("IterInversion")



        buffk =np.fromfile("Kr"+sufix+".bin",dtype='f4')
        Kr = np.reshape(buffk,[self.par.nodz,self.par.nody,self.par.nodx])
        Kr[-1,:,:] =  Kr[-2,:,:]

        
        buffk = np.fromfile("Kvp"+sufix+".bin",dtype='f4')
        Kvp = np.reshape(buffk,[self.par.nodz,self.par.nody,self.par.nodx])
        Kvp[-1,:,:] =  Kvp[-2,:,:]

        
        buffk = np.fromfile("Kvs"+sufix+".bin",dtype='f4')
        Kvs = np.reshape(buffk,[self.par.nodz,self.par.nody,self.par.nodx])
        Kvs[-1,:,:] =  Kvs[-2,:,:]
        
        subprocess.call("rm -f Kr" +sufix+".bin",shell=True)
        subprocess.call("rm -f Kvp" +sufix+".bin",shell=True)
        subprocess.call("rm -f Kvs" +sufix+".bin",shell=True)


        buffk = np.fromfile("PcondA"+sufix+".bin",dtype='f4')
        P = np.reshape(buffk,[self.par.nodz,self.par.nody,self.par.nodx])
        subprocess.call("rm -f PcondA"+sufix+".bin",shell=True)



        # Filtering Kernels
        Kr = filtK(Kr,D1,self.par)
        Kvp = filtK(Kvp,D1,self.par)
        Kvs = filtK(Kvs,D1,self.par)
        #P = filtK(P,D2,self.par)

        Kr = dampingK(Kr,2.0)
        Kvp = dampingK(Kvp,2.0)
        Kvs = dampingK(Kvs,2.0)

        Kr = filtK(Kr,D2,self.par)
        Kvp = filtK(Kvp,D2,self.par)
        Kvs = filtK(Kvs,D2,self.par)
        
        
        # Saving Files
        g = np.concatenate((Kr.flatten(),Kvp.flatten(),Kvs.flatten()))
        #g = np.concatenate((Kvp.flatten(),Kvs.flatten()))

        # if iter > 1 :
        #     g += e_m

        g[0:self.par.nodz * self.par.nody * self.par.nodx] = 0.0
            

        self.max_g = np.max(np.abs(g))
        
        g.astype('f4').tofile("g-"+sufix+".bin")
        sx = int(round(self.par.nodx/2)) 
        sy = int(round(self.par.nody/2))
        sz = int(round(self.par.nodz/2))
        
        plt.rc('font', size=5)
        Simshow(Kr,self.cx,self.cy,self.cz,sx,sy,sz,'Smooth_Krho')
        plt.savefig("Smth_Krho_it-"+str(iter)+".pdf",dpi=300,bbox_inches='tight')

        plt.rc('font', size=5)
        Simshow(Kvp,self.cx,self.cy,self.cz,sx,sy,sz,'Smooth_Kvp')
        plt.savefig("Smth_Kvp_it-"+str(iter)+".pdf",dpi=300,bbox_inches='tight')

        # plt.figure()
        # plt.imshow(D1)
        # plt.colorbar()
        # plt.savefig("filter.pdf",dpi=300,bbox_inches='tight')


        plt.rc('font', size=5)
        Simshow(Kvs,self.cx,self.cy,self.cz,sx,sy,sz,'Smooth_Kvs')
        plt.savefig("Smth_Kvs_it-"+str(iter)+".pdf",dpi=300,bbox_inches='tight')
        os.chdir(self.par.InvPar["-WorkDir"])



def filtK(K,D,par):
    nps=int(round((par.nodx * 0.05)))
    D1 = np.pad(D,((nps,nps),(nps,nps),(nps,nps)),'edge')
    K1 = np.pad(K,((nps,nps),(nps,nps),(nps,nps)),'edge')
    Kc = np.fft.fftshift(np.fft.fftn(K1)) * D1
    norm = 1.0 #par.nodx * par.nody * par.nodz
    K = np.real(np.fft.ifftn(np.fft.fftshift(Kc)))[nps:par.nodz+nps,nps:par.nody+nps,nps:par.nodx+nps]
    K *=norm
    return K



def KernelSmoothing(minlength,par):
    [zg ,yg, xg] = np.meshgrid(np.arange(0,par.nodz),np.arange(0,par.nody),\
                               np.arange(0,par.nodx),indexing='ij')
        
    f0x = 1.0 / (par.nodx * par.dx)
    f0y = 1.0 / (par.nody * par.dy)
    f0z = 1.0 / (par.nodz * par.dz)
    
    fc = 1 / minlength

    sig_x = int(round(fc / f0x))
    sig_y = int(round(fc / f0y))
    sig_z = int(round( fc / f0z))

    #print sig_x,sig_y,sig_z

    x0 = par.nodx/2
    y0 = par.nody/2
    z0 = par.nodz/2
    
    #Gaussian Function Kernel Smoothing                
    D1 = np.exp(-4.0 * (((xg-x0)**2.0) / (2 * sig_x * sig_x) +\
                        ((yg-y0)**2.0) / (2 * sig_y * sig_y) +\
                        ((zg-z0)**2.0) / (2 * sig_z * sig_z)))
    
    D1 = D1 / np.max(D1)
    
    return D1



def dampingK(kx,nstd):
    kx_mean = np.mean(kx)
    kx_std = np.std(kx)

    lwb = kx_mean - kx_std * nstd
    upb = kx_mean + kx_std * nstd
    
    kx[kx<lwb] = lwb 
    kx[kx>upb] = upb

    return kx
            





