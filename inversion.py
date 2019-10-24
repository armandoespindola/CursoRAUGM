import os,time
import subprocess
import numpy as np
from inversion.inversion import inversion


TimeSimulation=4
d_freq = 2.0
IterPerFreq = 11
FreqMax = 3.0
nf=2
TypeAdj=1
veldisp = 'disp'
slurm = 0
# Data Format
syn_sac = False
obs_sac = False

# Smoothing mlen1 -- Kernel  mlen2 --> Precond

mlen1 = 1000.0  # [m]
mlen2 = 1000.0 # [m]

freq = np.arange(1.0,FreqMax + d_freq,d_freq)
Inv = inversion(IterPerFreq,TypeAdj,veldisp,slurm,mlen1,mlen2,syn_sac,obs_sac)

print freq

path = os.getcwd()

# Erase Iterations Directory 
subprocess.call("rm -r IterInversion/*",shell=True)

kf = 0
for j in range(0,freq.size-nf+1,nf-1):
    #print 'freq:',j
    kf +=1

    # CREATES FREQUENCY FILE 
    fstat  = open('frequency.par','w')
    fstat.write('// FREQUENCIES\n')
    
    for fi in range(0,nf):
        #print freq[fi+j]
        fstat.write('%f\n' % (freq[fi+j]))

    fstat.close()

    # # STARTS INVERSION

    Inv.Start()

    os.chdir(path)

    subprocess.call("mkdir -p IterInversion/f"+str(kf),shell=True)
    subprocess.call("cp IterInversion/m*.bin IterInversion/f"+str(kf) +"/",shell=True)
    subprocess.call("cp IterInversion/*.pdf IterInversion/f"+str(kf) + "/",shell=True)
    subprocess.call("cp IterInversion/*.dat IterInversion/f"+str(kf) + "/",shell=True)
