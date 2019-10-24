import numpy as np
import subprocess
import os
import shutil
import time
import sys



class solver:
        def __init__(self,par):
                self.par = par

	def Compile(self):
		os.chdir(self.par.InvPar["-WorkDir"]);
		os.chdir("wave3Dfd/");
		subprocess.call(["make", "clean"]);
		subprocess.call(["make", "all"]);

    
	def CopySolver(self):
		os.chdir(self.par.InvPar["-WorkDir"]);
		os.chdir("wave3Dfd/");
		for i in range(0,self.par.EventNum):
			dest = self.par.InvPar["-WorkDir"]+"/"+self.par.EventName[i]+ \
			"/"+self.par.EventName[i]+".out";
			shutil.copy("bin/wave3Dfd.out",dest);

	def Forward(self):
		os.chdir(self.par.InvPar["-WorkDir"]);
		for i in range(0,self.par.EventNum):
			program = "mpirun -np " + self.par.InvPar["-nmpi"] + \
			 " ./" + self.par.EventName[i] +".out" \
			+ " -nFile " +self.par.InvPar["-parForw"] + \
			" > ForwardModeling.log";
			print "Running (Forward Simulation): " + self.par.EventName[i] +".out" \
			+ " -nFile " +self.par.InvPar["-parForw"];
			os.chdir(self.par.EventName[i]);
			subprocess.call(program,shell=True);
			os.chdir("../")

	def Adjoint(self):
		os.chdir(self.par.InvPar["-WorkDir"]);
		for i in range(0,self.par.EventNum):
			program = "mpirun -np " + self.par.InvPar["-nmpi"] + \
			 " ./" + self.par.EventName[i] +".out" \
			+ " -nFile " +self.par.InvPar["-parAdj"] + \
			" > AdjointModeling.log";
			print "Running (Adjoint Simulation) : " +  self.par.EventName[i] +".out" \
			+ " -nFile " +self.par.InvPar["-parAdj"];
			os.chdir(self.par.EventName[i]);
			subprocess.call(program,shell=True);
			os.chdir("../")


        def job_status(self,job):
                """ Queries completion status of a single job
                """
                stdout = ""
                while not stdout:
                        stdout = subprocess.check_output('sacct -n -o jobid,state -j '+ job.split('_')[0],
                                                         shell=True)

                state = "WAITING"
                lines = stdout.strip().split('\n')
                #print lines,'line'
                for line in lines:
                        if line.split()[0]==job:
                                state = line.split()[1]
				#print job,'\t',state
                return state

        def job_array_status(self,jobs):
                """ Determines completion status of job or job array
                """
                states = []

                for job in jobs:
                        state = self.job_status(job)
                        if state in ['TIMEOUT']:
                                print "TIMEOUT"
                                sys.exit(-1)
                        elif state in ['FAILED', 'NODE_FAIL']:
                                print "FAILED"
                                sys.exit(-1)
                        elif state in ['COMPLETED']:
                                states += [1]
                        else:
                                states += [0]

                isdone = all(states)

                return isdone, jobs


        def job_id_list(self,stdout, ntask):
                """ Parses job id list from sbatch standard output
                """
                job_id = stdout.split()[-1].strip()

		num =[]

		for jj in range(0,self.par.EventNum):
			#print jj
			num += [''.join([ii for ii in self.par.EventName[jj] if ii.isdigit()])]

                return [job_id+'_'+ kk for kk in num]




        def SlurmForward(self):
                os.chdir(self.par.InvPar["-WorkDir"])
                nshot = self.par.EventNum

                stdout = subprocess.check_output('sbatch Forward.sh',shell=True)
                print stdout

                jobs = self.job_id_list(stdout,nshot)

                while True:
                        # wait a few seconds between queries
                        time.sleep(5)
			b = "RUNNING FORWARD" + "...\r" 
			print b

                        isdone, jobs = self.job_array_status(jobs)
                        #print isdone
                        if isdone:
                                print "Forward Done"
                                break







        def SlurmAdjoint(self):
                os.chdir(self.par.InvPar["-WorkDir"])
                nshot = self.par.EventNum

                stdout = subprocess.check_output('sbatch Adjoint.sh',shell=True)
                print stdout

                jobs = self.job_id_list(stdout,nshot)

                while True:
                        # wait a few seconds between queries
                        time.sleep(5)
			b = "RUNNING ADJOINT" + "...\r"
                        print b
                        isdone, jobs = self.job_array_status(jobs)
                        #print isdone
                        if isdone:
                                print "Adjoint Done"
                                break


	



        
                
                
                




