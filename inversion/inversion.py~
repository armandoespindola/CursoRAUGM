import os
import subprocess

from Solver.solver import solver
from Tools.parameter import param
from Tools.UpdateModel import UpdateParam
from Postprocess.ProcKernel import ProcKernel
from AdjointSource.AdjointSourceHarmonic import AdjointSource
from ForwardSource.ForwardSource import ForwardSource
from Optimization.optimization import optimization as optim
import time,sys
import numpy as np


class inversion:
    def __init__(self,IterFreq,TypeAdj,veldisp,slurm,mlen1,mlen2,syn_sac,obs_sac):
        # Smoothing Size
        # Kernel --> mlen1
        # Precond --> mlen2
        self.iter = IterFreq
        self.ntimes = 0
        self.TypeAdj = TypeAdj
        self.veldisp = veldisp
        self.slurm = slurm
        self.mlen1 = mlen1
        self.mlen2 = mlen2
        self.syn_sac = syn_sac
        self.obs_sac = obs_sac

    def Start(self):

        self.ntimes += 1

        # Parameters Simulation 
        par = param()

        #### Initialize Work Envioroment
        ####
        ####

        par.InitEnviron()
        E = AdjointSource(par,self.veldisp,self.syn_sac,self.obs_sac)
        FwdSrc = ForwardSource(par)
        opt = optim()
        ModelPar = UpdateParam(par)

        if (self.ntimes > 1):
            ### Return back to new model
            ModelPar.NewRestart()

        # Simulation Initial Iteration
        solv = solver(par)

        if (self.ntimes <= 1):
            solv.Compile()
            solv.CopySolver()

        ### SOLVE FORWARD MODELING

        if self.slurm == 0:
            solv.Forward()
        elif self.slurm == 1 :
            solv.SlurmForward()
            
        #sys.exit(0)

        ### COMPUTE ADJOINT SOURCE ADN L2NORM
        E.AdjSrc(self.TypeAdj,0)
        E.Error(self.TypeAdj,"new",0)
        #sys.exit(0)


        # Copy Initial Model as New Only First Iteration
        ModelPar.UpdateNew()

        # Optimization Pamameters
        K = ProcKernel(par)


        # Initialization Optimization Method
        opt.Setup(par)


        for iterNCG in range(1,self.iter+1):
            print "##############"
            print " Iter: ", iterNCG
            print "##############"

            ### SOLVE ADJOINT
            FwdSrc.FwdSrc()
            if self.slurm == 0:
                solv.Adjoint()
            elif self.slurm == 1:
                solv.SlurmAdjoint()
            #sys.exit(0)

            # Kernel Processing
            K.SumKernel("new",iterNCG)
            K.Ksmoth(iterNCG,"new",self.mlen1,self.mlen2)

            #sys.exit(0)         
            
            print "MAX_VALUE_GRAD: ",K.max_g

            #if K.max_g < 5.0:
            #    print "Small Changes in Gradient"
            #    break
            # sys.exit(0)
        
            ### COMPUTE GRADIENT
            opt.ComputeDirection()

            ### LINE SEARCH
            linesearch(ModelPar,opt,solv,E,self.TypeAdj,self.slurm,iterNCG)

            ### SOLVE FORWARD MODELING
            #solv.SlurmForward()
            #solv.Forward()
        
            ##### SOLVE ADJOINT
            E.AdjSrc(self.TypeAdj,iterNCG)
            E.Error(self.TypeAdj,"new",iterNCG)





def linesearch(ModelPar,opt,solv,E,TypeAdj,slurm,iter):
    
    ### INITIALIZE LINE SERARCH
    opt.InitSearch()

    #### LINE SEARCH
    while True:
        #### UPDATE FILES VP.bin RHO.bin VS.bin WITH M.TEST
        error =  ModelPar.UpdateTest()

        if error==0:
            print "Test-Model Update"
        #### SOLVE FORWARD PROBLEM
            if slurm == 0:
                solv.Forward()
            elif slurm == 1:
                solv.SlurmForward()

            #### COMPUTE ERROR (MISFIT)
            E.Error(TypeAdj,"test",iter)

        else:
            ferror = np.loadtxt(E.par.FoldIter+'/f-new.txt') * 1.1
            f = open(E.par.FoldIter+'/f-'+'test'+'.txt', 'w')
            f.write('%e\n' % ferror)
            f.close()
        #sys.exit(0)

        #print " trial step", opt.opt.LineSearch.step_count + 1
        status = opt.UpdateSearch()

        if status == 1:
            error = ModelPar.UpdateTest()
            opt.FinelizeSearch()
            break
        elif status < 0:
            if (opt.retry_status()):
                opt.restart()
                linesearch(ModelPar,opt,solv,E,TypeAdj,slurm,iter)
                break
            else:
                print ' Line search failed\n\n'
                sys.exit(-1)
        else:
            continue



