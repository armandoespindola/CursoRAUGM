import subprocess
import numpy as np
#from inversion.Tools.parameter import param

class UpdateParam:
	def __init__(self,par):
		self.par = par

	def UpdateTest(self):
                mtest = np.fromfile(self.par.FoldIter +"/m-test.bin",dtype='f4')
                mnew = np.fromfile(self.par.FoldIter +"/m-new.bin",dtype='f4')
                #rhon = np.fromfile(self.par.FoldIter +"/rho-new.bin",dtype='f4')
                nstride = self.par.nodx * self.par.nody * self.par.nodz

                rho = mtest[0:nstride]
                vp = mtest[nstride:2 * nstride]
                vs = mtest[2 * nstride:3 * nstride]

                #vp = mtest[0:nstride]
                #vs = mtest[nstride:2 * nstride]

                #vpn = mnew[0:nstride]
                #vsn = mnew[nstride:2 * nstride]

                #rho = rhon + (mtest[nstride:2 * nstride] - mnew[nstride:2 * nstride]) * 0.0

                ModelOK = self.TestModel(vp,vs,rho)
                
                #vp = mtest[nstride:2 * nstride]
                #vs = mtest[2 * nstride:3 * nstride]

                

                #rhon = mnew[0:nstride]
                #vpn = mnew[nstride:2 * nstride]
                #vsn = mnew[2 * nstride:3 * nstride]

                # a1 = vp/vs < 1.6
		# a2 = vp/vs > 1.9
                
		# rho[a1] = rhon[a1]
		# vp[a1] = vpn[a1]
               	# vs[a1] = vsn[a1]
		
		# rho[a2] = rhon[a2]
                # vp[a2] = vpn[a2]
                # vs[a2] = vsn[a2]

                #rho = rhon
                #vp = vpn
                #vs = vsn
                
                if ModelOK==0:
                        print "Good Model"
                        rho.astype('f4').tofile(self.par.WorkDir+"/RHO.bin") 
                        vp.astype('f4').tofile(self.par.WorkDir+"/VP.bin")
                        vs.astype('f4').tofile(self.par.WorkDir+"/VS.bin")
                else:
                        print "Error Model"
                return ModelOK

        def NewRestart(self):
                mtest = np.fromfile(self.par.FoldIter +"/m-new.bin",dtype='f4')
                #rho  = np.fromfile(self.par.FoldIter +"/rho-new.bin",dtype='f4')
                nstride = self.par.nodx * self.par.nody * self.par.nodz

                #vp = mtest[0:nstride]
                #vs = mtest[nstride:2 * nstride]

                rho = mtest[0:nstride]
                vp = mtest[nstride:2 * nstride]
                vs = mtest[2 * nstride:3 * nstride]

                rho.astype('f4').tofile(self.par.WorkDir+"/RHO.bin")
                vp.astype('f4').tofile(self.par.WorkDir+"/VP.bin")
                vs.astype('f4').tofile(self.par.WorkDir+"/VS.bin")
                

	def UpdateNew(self):
                rho = np.fromfile(self.par.WorkDir +"/RHO.bin",dtype='f4')
                vp = np.fromfile(self.par.WorkDir +"/VP.bin",dtype='f4')
                vs = np.fromfile(self.par.WorkDir +"/VS.bin",dtype='f4')
                m = np.concatenate((rho,vp,vs))
                #m = np.concatenate((vp,vs))
                m.astype('f4').tofile(self.par.FoldIter +"/m-new.bin")
                #rho.astype('f4').tofile(self.par.FoldIter +"/rho-new.bin")

        def TestModel(self,vp,vs,rho):
                mu = rho * vs**2
                lamb = (rho * vp**2) - 2.0 * mu
                error = 0
                C0 = 1.0 / 24.0;
                C1 = 9.0 / 8.0;

                if (np.any(mu < 0)) or (np.any(rho <= 0)) or (np.any(lamb <=0)) :
                        error = 1
                        print "ERROR MODEL PARAMETERS"
                        return error

                kx = self.par.dt * np.sqrt(3.0) * vp *(C0 + C1) / self.par.dx
                ky = self.par.dt * np.sqrt(3.0) * vp *(C0 + C1) / self.par.dy
                kz = self.par.dt * np.sqrt(3.0) * vp *(C0 + C1) / self.par.dz
                
                if (np.any(kx >= 1)) or (np.any(ky >= 1)) or (np.any(kz >= 1)):
                        error = 1
                        print "CFL NOT SATISFIED"
                        return error

                # if (np.any(vp/vs < 1.55)):
                #         error = 1
                #         print "VP/VS ratio out of bounds"
                #         return error
                        
                return error
                
