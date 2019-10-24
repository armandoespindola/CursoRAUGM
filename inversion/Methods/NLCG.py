import os
import numpy as np
import subprocess
from inversion.Tools.parameter import param
from inversion.Tools.UtilFunc import dot,save

class NLCG(param):
	""" Nonlinear conjugate gradient method
	"""
	def __init__(self,maxiter = np.inf,thresh = 1.):

		param.__init__(self);
		self.Fiter = self.FoldIter+'/iterNCLG';
		os.chdir(self.FoldIter);
		subprocess.call(["mkdir", "-p",self.Fiter]);
		self.maxiter = maxiter;
		self.thresh = thresh;
		self.iter = 0;
		os.chdir(self.WorkDir);

	def __call__(self):
		os.chdir(self.FoldIter)
		
		self.iter += 1;
		save(self.FoldIter+'/'+'iterNCLG/iter.txt', self.iter);

		gnew = np.fromfile("g-new.bin",dtype='f4')

		if self.iter == 1:
			return -gnew, 0

		elif self.iter > self.maxiter:
			print 'restarting NLCG... [periodic restart]'
			self.restart()
			return -gnew, 1   
           # compute search direction
		gold = np.fromfile("g-old.bin",dtype='f4')
		pold = np.fromfile("p-old.bin",dtype='f4')
  
		beta = pollak_ribere(gnew, gold)
		pnew = -gnew + beta*pold


        # check restart conditions
		if check_conjugacy(gnew, gold) > self.thresh:
                  print 'restarting NLCG... [loss of conjugacy]'
                  self.restart()
                  return -gnew, 1

		elif check_descent(pnew, gnew) > 0.:
                  print 'restarting NLCG... [not a descent direction]'
                  self.restart()
                  return -gnew, 1
		else:
                  return pnew, 0
        
	def restart(self):
         """ Restarts algorithm
         """
         self.iter = 1
         save(self.FoldIter+'/'+'iterNCLG/iter.txt', self.iter)



def pollak_ribere(gnew, gold):
        num = dot(gnew, gnew-gold)
        den = dot(gold, gold)
        beta = num/den
        return beta

def check_conjugacy(gnew, gold):
        a = abs(dot(gnew, gold) / dot(gnew, gnew))
        print "check_conjugancy:",a
        return a

def check_descent(pnew, gnew):
        a = dot(pnew, gnew) / dot(gnew, gnew)
        print "check_decent:",a
        return a














