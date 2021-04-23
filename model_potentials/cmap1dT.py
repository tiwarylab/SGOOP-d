# -------------------------------------------------------------------------------------------------
# The class is used for estimating SGOOP-based kinetic distances.
# -------------------------------------------------------------------------------------------------

import numpy as np

class SgoopDistance:
    
    def __init__(self, coefficient, eVal, eVec, binedges):
        
        self.coefficient = coefficient
        self.eVal = eVal
        self.eVec = eVec
        self.binedges = binedges
        
        self.rc = {}
        self.idx = {}
        
    
    def proj2rc(self, state_pos):
        """
        state_pos: (x1, x2, x3, ..., xN).
        coefficient: (c1, c2, c3, ... cN).
        This will return rc value with rc=c1*x1+c2*x2+c3*x3+...+cN*xN.
        """
    
        if len(state_pos)!=len(self.coefficient):
            raise Exception('state position vector should have the same length as coefficient of RC.')
    
        return np.dot(self.coefficient, state_pos)

    def pairwise_d(self, pos1, pos2, N=None):
        """
        pos1, pos2: The state position values of two states. pos_i is ususally a set of colvars.
        eVal, eVec: The eigenvalues and eigenvectors of maximum caliber-based rate matrix.
        binedges: The binedges of SGOOP.
        """
        
        rc1, rc2 = self.proj2rc(pos1), self.proj2rc(pos2)
        self.rc['1'] = rc1
        self.rc['2'] = rc2
        
        if not isinstance(rc1, float) or not isinstance(rc2, float): 
            raise Exception('rc1 or rc2 is not float or numpy.float64')
        if N is None:
            N=self.eVal.shape[0]
        self.binedges[-1] += 1000 # Add 1000 to ensure that the right edge is inclusive.
    
        idx1, idx2 = np.digitize([rc1, rc2], self.binedges)-1
        self.idx['1'] = idx1
        self.idx['2'] = idx2
        
        rate = -self.eVal
        d_comm = [(0.5/rate[k])*(self.eVec[:,k][idx1]/self.eVec[:,0][idx1]\
                                      -self.eVec[:,k][idx2]/self.eVec[:,0][idx2])**2 for k in range(1,N)]
        d_comm = np.sum(d_comm)
    
        return d_comm
    