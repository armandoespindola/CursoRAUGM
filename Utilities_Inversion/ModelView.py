import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def colormy():
    colors = [(0.50,0,0),(1,0,0),(1.0,1.0,0.0),(1,1,1),(0,1,1),(0,0,1),(0,0,0.56)] 
    #colors = colors[::-1]
    cmap_name = "fwi"
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=21)
    return cm


nx = 128
ny = 128
nz = 64

buff =np.fromfile("TargetModel/VP.bin",dtype='f4')
vp_t = np.reshape(buff,[nz,ny,nx])

vp_min = np.min(vp_t)
vp_max = np.max(vp_t)


buff =np.fromfile("InitialModel/VP.bin",dtype='f4')
vp0 = np.reshape(buff,[nz,ny,nx])

buff =np.fromfile("VP.bin",dtype='f4')
vp_inv = np.log(np.reshape(buff,[nz,ny,nx])/vp0)


buff =np.fromfile("TargetModel/VS.bin",dtype='f4')
vs_t = np.reshape(buff,[nz,ny,nx])

vs_min = np.min(vs_t)
vs_max = np.max(vs_t)

buff =np.fromfile("InitialModel/VS.bin",dtype='f4')
vs0 = np.reshape(buff,[nz,ny,nx])

buff =np.fromfile("VS.bin",dtype='f4')
vs_inv = np.log(np.reshape(buff,[nz,ny,nx])/vs0)

cmap = colormy()


plt.figure()
plt.subplot(3,1,1)
plt.pcolormesh(vp0[:,:,64],cmap=cmap,vmin = vp_min,vmax = vp_max)
plt.colorbar()
plt.subplot(3,1,2)
plt.pcolormesh(vp_inv[:,:,64],cmap=cmap,vmin = np.min(vp_inv),vmax = -np.min(vp_inv))
plt.colorbar()
plt.subplot(3,1,3)
plt.pcolormesh(vp_t[:,:,64],cmap=cmap,vmin = vp_min,vmax = vp_max)
plt.colorbar()

plt.figure()
plt.subplot(3,1,1)
plt.pcolormesh(vs0[:,:,64],cmap=cmap,vmin = vs_min,vmax = vs_max)
plt.colorbar()
plt.subplot(3,1,2)
plt.pcolormesh(vs_inv[:,:,64],cmap=cmap,vmin = -np.max(vs_inv),vmax = np.max(vs_inv))
plt.colorbar()
plt.subplot(3,1,3)
plt.pcolormesh(vs_t[:,:,64],cmap=cmap,vmin = vs_min,vmax = vs_max)
plt.colorbar()


plt.figure()
plt.subplot(3,1,1)
plt.pcolormesh(vp0[:,64,:],cmap=cmap,vmin = vp_min,vmax = vp_max)
plt.colorbar()
plt.subplot(3,1,2)
plt.pcolormesh(vp_inv[:,64,:],cmap=cmap,vmin = np.min(vp_inv),vmax = -np.min(vp_inv))
plt.colorbar()
plt.subplot(3,1,3)
plt.pcolormesh(vp_t[:,64,:],cmap=cmap,vmin = vp_min,vmax = vp_max)
plt.colorbar()

plt.figure()
plt.subplot(3,1,1)
plt.pcolormesh(vs0[:,64,:],cmap=cmap,vmin = vs_min,vmax = vs_max)
plt.colorbar()
plt.subplot(3,1,2)
plt.pcolormesh(vs_inv[:,64,:],cmap=cmap,vmin = -np.max(vs_inv),vmax = np.max(vs_inv))
plt.colorbar()
plt.subplot(3,1,3)
plt.pcolormesh(vs_t[:,64,:],cmap=cmap,vmin = vs_min,vmax = vs_max)
plt.colorbar()

plt.show()
