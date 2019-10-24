# Program Calls ForwardSolvers For all Events
# Wed Apr 17 23:34:30 +03 2019

import os,sys
import subprocess
sys.path.append('inversion')
print os.getcwd()
from Solver.solver import solver
from Tools.parameter import param
from Tools.UpdateModel import UpdateParam



# Parameters Simulation 
par = param()

ModelPar = UpdateParam(par)

### Return back to new model
ModelPar.NewRestart()

