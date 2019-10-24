clear
clc
close all

 fidx=fopen('../InitialModel/RHO.bin','r');
 rho=fread(fidx,'float');
 fidx=fopen('../InitialModel/VP.bin','r');
 vp=fread(fidx,'float');
 fidx=fopen('../InitialModel/VS.bin','r');
 vs=fread(fidx,'float');
 

 fidx=fopen('../RHO.bin','r');
 rho_inv=fread(fidx,'float')- rho;
 fidy=fopen('../TargetModel/RHO.bin','r');
 rho0=fread(fidy,'float');%- rho;
 
 %fidx=fopen('../../../VP.bin','r');
 fidx=fopen('../IterInversion/VP.bin','r');
 vp_inv=fread(fidx,'float');%-vp;
 fidy=fopen('../TargetModel/VP.bin','r');
 vp0=fread(fidy,'float');%- vp;
 
 %fidx=fopen('../../../VS.bin','r');
 fidx=fopen('../IterInversion/VS.bin','r');
 vs_inv=fread(fidx,'float');% - vs;
 fidy=fopen('../TargetModel/VS.bin','r');
 vs0=fread(fidy,'float');%-vs;
 
ni=128;
nj=128;
nk=64;

Mrho_inv= permute(reshape(rho_inv,ni,nj,nk),[2,1,3]);
Mrho0= permute(reshape(rho0,ni,nj,nk),[2,1,3]);
Mvp_inv= permute(reshape(vp_inv,ni,nj,nk),[2,1,3]);
Mvp0= permute(reshape(vp0,ni,nj,nk),[2,1,3]);
Mvs_inv= permute(reshape(vs_inv,ni,nj,nk),[2,1,3]);
Mvs0= permute(reshape(vs0,ni,nj,nk),[2,1,3]);

figure()
slice(Mrho0,64,64,32);
title('RHO TARGET')
colormap(jet)
colorbar
shading interp
figure()
slice(Mrho_inv,64,64,32);
title('RHO INV')
colormap(jet)
colorbar
shading interp

figure()
slice(Mvp0,64,64,32);
colormap(redblue)
title('VP TARGET')
colormap(jet)
colorbar
shading interp
figure()
slice(Mvp_inv,64,64,32);
title('VP INV')
colormap(jet)
colorbar
shading interp

figure()
slice(Mvs0,64,64,32);
title('VS TARGET')
colormap(jet)
colorbar
shading interp
figure()
slice(Mvs_inv,64,64,32);
title('VS INV')
colormap(jet)
colorbar
shading interp
