import numpy as np
import os
import subprocess

class param:
	def __init__(self):
		
		# Event Folders
		self.EventName = np.genfromtxt('Events.par',dtype='str',
			comments='#',delimiter=' ');
		
		self.EventNum = self.EventName.size;

		self.InvPar = {};
		
		# Inversion Parameters
		par = np.genfromtxt('InvParam.par',dtype='str',
			comments='#',delimiter=' ');

		for i in range(0,par.shape[0]):
			self.InvPar.update({par[i][0]:par[i][1]})

                self.InvPar.update({"-WorkDir":os.getcwd()})

		self.WorkDir = self.InvPar["-WorkDir"];

		self.FoldIter = self.InvPar["-WorkDir"] + '/IterInversion';

		self.ModPar = {};
		
		# Modeling Parameters
		mpar = np.genfromtxt('ParForward.par',dtype='str',
			comments='//',delimiter=' ');
		for i in range(0,mpar.shape[0]):
			self.ModPar.update({mpar[i][0]:mpar[i][1]})
		self.xi = float(self.ModPar["-xi"]);
		self.xf = float(self.ModPar["-xf"]);
		self.yi = float(self.ModPar["-yi"]);
		self.yf = float(self.ModPar["-yf"]);
		self.zi = float(self.ModPar["-zi"]);
		self.zf = float(self.ModPar["-zf"]);
		self.nodx = int(self.ModPar["-nx"]);
		self.nody = int(self.ModPar["-ny"]);
		self.nodz = int(self.ModPar["-nz"]);
		self.f0 = float(self.ModPar["-f0"]);
		self.nsrc = int(self.ModPar["-ns"]);
		self.nrec = int(self.ModPar["-nr"]);
		self.dt = float(self.ModPar["-dt"]);
		self.Tsim = float(self.ModPar["-t"]);
                self.nt = int(self.Tsim / self.dt)
		self.dx = (self.xf-self.xi) / (self.nodx - 1);
		self.dy = (self.yf-self.yi) / (self.nody - 1);
		self.dz = (self.zf-self.zi) / (self.nodz - 1);

                self.freq = np.genfromtxt('frequency.par',dtype='f4',
                                          comments='//',delimiter=' ');

	def InitEnviron(self):
		subprocess.call(["mkdir", "-p",self.FoldIter]);
		os.chdir(self.FoldIter);
		#subprocess.call("rm -f *.*",shell=True);
		os.chdir(self.InvPar["-WorkDir"]);
