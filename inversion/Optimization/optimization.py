import numpy as np 
from inversion.Optimization.NLCG import NLCG
from inversion.Tools.parameter import param
from inversion.Tools.UtilFunc import polyfitB,angle,dot,backtrack2,count_zeros,save
import subprocess
import os,sys

import math


class optimization:
	def Setup(self,par):
		self.par = par
                os.chdir(self.par.InvPar["-WorkDir"])
		self.maxiter = int(self.par.InvPar["-maxiter"]) 
		self.thresh = float(self.par.InvPar["-thresh"])
		self.step_count_max = int(self.par.InvPar["-stepcountmax"])
		self.step_len_max = float(self.par.InvPar["-steplenmax"])
		self.stepInit = int(self.par.InvPar["-StepMaxAuto"])
		
		# Optimization Method Initialization

                if (self.par.InvPar["-method"] == "LBFGS"):
                        print " Optimization Method LBFGS "
                        print " Not implemented"
                        # self.opt = LBFGS();
                        # self.opt.Setup(self.maxiter,self.thresh)
                        # self.opt.Setup_LineSearch(self.step_count_max,self.step_len_max)
                        
                elif (self.par.InvPar["-method"] == "NLCG"):
                        print " Optimization Method NLCG "
                        self.opt = NLCG();
                        self.opt.Setup(self.maxiter,self.thresh)
                        self.opt.Setup_LineSearch(self.step_count_max,self.step_len_max)
                else:
                        print " No optimization Method Selected "
                        sys.exit(0)


# Compute Drirection
	def ComputeDirection(self):
		self.opt.Compute_Direction()



	def InitSearch(self):
		""" Computes step length in line search
		"""

		
		m = np.fromfile(self.par.FoldIter+'/m-new.bin',dtype='f4')
		g = np.fromfile(self.par.FoldIter+'/g-new.bin',dtype='f4')
		p = np.fromfile(self.par.FoldIter+'/p-new.bin',dtype='f4')
		f = np.loadtxt(self.par.FoldIter+'/f-new.txt')
		norm_m = np.max(np.abs(m))
		norm_p = np.max(np.abs(p))
		gtg = dot(g,g)
		gtp = dot(g,p)

		if self.opt.restarted:
			self.opt.LineSearch.clear_history()


		#deltaL,_ = self.opt.LineSearch.initialize(0,f,gtg,gtp)


                print "LEN_MAX1:",self.opt.LineSearch.step_len_max

		if self.stepInit:
			self.opt.LineSearch.step_len_max = 0.05 * norm_m / norm_p

                print "LEN_MAX2:",self.opt.LineSearch.step_len_max

                deltaL,_ = self.opt.LineSearch.initialize(0,f,gtg,gtp)

	
		if self.stepInit and len(self.opt.LineSearch.step_lens)<=1:
			deltaL = self.step_len_max * self.opt.LineSearch.step_len_max


		save(self.par.FoldIter+'/deltaL',deltaL)

		m_test = m + deltaL * p;

		m_test.astype('f4').tofile(self.par.FoldIter+'/m-test.bin')


	def UpdateSearch(self):
		""" Updates line search step length

		Program status:
			status > 0 : finished
			status == 0 : not finished
			status <  0 : fail		
		"""

		delt = np.loadtxt(self.par.FoldIter+'/deltaL');
		ftest = np.loadtxt(self.par.FoldIter+'/f-test.txt')


		deltaL,status = self.opt.LineSearch.update(delt,ftest)

		#print "status:",status,deltaL,self.opt.LineSearch.step_len_max

		if status >= 0:
			m = np.fromfile(self.par.FoldIter+'/m-new.bin',dtype='f4')
			p = np.fromfile(self.par.FoldIter+'/p-new.bin',dtype='f4')

			save(self.par.FoldIter+'/deltaL',deltaL)
			
			m_test = m + deltaL * p;
			m_test.astype('f4').tofile(self.par.FoldIter+'/m-test.bin')

		return status



	def FinelizeSearch(self):

		m = np.fromfile(self.par.FoldIter+'/m-new.bin',dtype='f4')
		g = np.fromfile(self.par.FoldIter+'/g-new.bin',dtype='f4')
		p = np.fromfile(self.par.FoldIter+'/p-new.bin',dtype='f4')
		f = self.opt.LineSearch.search_history()[1]
		deltas = self.opt.LineSearch.search_history()[0]


		if self.opt.iter() > 1:
			command = "rm "+ self.par.FoldIter+'/m-old.bin'
			subprocess.call(command,shell=True);

			command = "rm "+ self.par.FoldIter+'/g-old.bin'
			subprocess.call(command,shell=True);

			command = "rm "+ self.par.FoldIter+'/p-old.bin'
			subprocess.call(command,shell=True);

			command = "rm "+ self.par.FoldIter+'/f-old.txt'
			subprocess.call(command,shell=True);


		os.chdir(self.par.FoldIter)

                command = "mv m-new.bin m-old.bin"
		subprocess.call(command,shell=True)

                command = "mv f-new.txt f-old.txt"
		subprocess.call(command,shell=True)

                command = "mv g-new.bin g-old.bin"
		subprocess.call(command,shell=True)

                command = "mv p-new.bin p-old.bin"
		subprocess.call(command,shell=True);

		command = "cp m-test.bin m-new.bin"
		subprocess.call(command,shell=True)
		save(self.par.FoldIter+'/f-new.txt',f.min())

		os.chdir(self.par.WorkDir)

		# Log Inversion

                fstat  = open(self.par.FoldIter+'/stats-'+str(self.opt.iter())+'.txt','w')

                fstat.write('factor : %e\n' % (-dot(g,g)**-0.5 * (f[1]-f[0])/(deltas[1]-deltas[0])))
                fstat.write('gradient_norm_L1: %e\n' % (np.linalg.norm(g, 1)))
                fstat.write('gradient_norm_L2: %e\n' % (np.linalg.norm(g, 2)))
                fstat.write('misfit: %e\n' % (f[0]))
                #fstat.write('restarted: ', self.opt.restarted)
                fstat.write('slope: %e\n' % ((f[1]-f[0])/(deltas[1]-deltas[0])))
                fstat.write('step_count: %d\n' %  (self.opt.LineSearch.step_count))
                fstat.write('step_length: %e\n' % (deltas[f.argmin()]))
                fstat.write('theta: %e\n' % (180.*np.pi**-1*angle(p,-g)))

                fstat.close()

                

        def restart(self):
                g = (-1.0) * np.fromfile(self.par.FoldIter+'/g-new.bin',dtype='f4')
                g.astype('f4').tofile(self.par.FoldIter+'/p-new.bin')
                self.opt.LineSearch.clear_history()
                self.opt.restarted = 1
                self.opt.LineSearch.writer.iter -= 1
                self.opt.LineSearch.writer.newline()

        def retry_status(self):
                print "Steplength - Retry"
                g = np.fromfile(self.par.FoldIter+'/g-new.bin',dtype='f4')
                p = np.fromfile(self.par.FoldIter+'/p-new.bin',dtype='f4')
                theta = angle(p,-g)
                thresh = 1.e-3
                
                if abs(theta) < thresh:
                        return 0
                else:
                        return 1





















		
