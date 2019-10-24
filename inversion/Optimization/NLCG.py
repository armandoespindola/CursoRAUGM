import numpy as np 
from inversion.Methods.NLCG import NLCG as m_NLCG
from inversion.LineSearch.bracketing import Bracket

class NLCG:
         """ Nonlinear conjugate gradient method
         """
         def Setup(self,maxiter,thresh):
                  self.CG = m_NLCG(maxiter,thresh)

         def Setup_LineSearch(self,step_count_max,step_len_max):
                  self.LineSearch = Bracket(step_count_max,step_len_max)


         def Compute_Direction(self):
                  pnew,self.restarted = self.CG()
                  pnew.astype('f4').tofile("p-new.bin");

         def Restart(self):
                  self.CG.restart()

         def iter(self):
                  return self.CG.iter