# Program Calls ForwardSolvers For all Events
# Wed Apr 17 23:34:30 +03 2019

import os,sys
import subprocess
sys.path.append('inversion')
print os.getcwd()
from Solver.solver import solver
from Tools.parameter import param



# Parameters Simulation 
par = param()

# Simulation Initial Iteration
solv = solver(par)
solv.Compile()
solv.CopySolver()

### SOLVE FORWARD MODELING
solv.Forward()
