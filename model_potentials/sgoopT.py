# -------------------------------------------------------------------------------------------------
# The class is used to perform spectral gap calculation from either biased or unbiased simulation.
# -------------------------------------------------------------------------------------------------

import numpy as np

def stationary(rc, rc_bin, weights=None):
    """
    Build stationary probability distribution.
    
    Input arguments are rc, rc_bin, weights (default None)
    rc: The trajectory of reaction coordinate value, 
    rc_bin: The bin number used for spatially discretizing rc.
    weights: If True, this will be used to reweight probability.
    
    The output is the probability vector.
    """
    if weights is None:
        hist, binedges=np.histogram(rc, bins=rc_bin)
        prob=hist.T/np.sum(hist.T)
    else:
        hist, binedges=np.histogram(rc, weights=weights, bins=rc_bin)
        prob=hist.T/np.sum(hist.T)
    
    return prob, binedges
    
def mu_factor(unbiased_rc, pi, rc_bin, max_d=1):
    """
    The dynamical prefactor in maxcal estimated transition rate matrix. 
    This was related to the constraints of the number of first nearest neighbor transitions.
    """
    rc_min, rc_max=np.min(unbiased_rc), np.max(unbiased_rc)
    bins=np.linspace(rc_min, rc_max, rc_bin+1)
    unbiased_rc_idx=np.digitize(unbiased_rc, bins)
    
    d=np.arange(max_d)+1
    MU=np.zeros(max_d)
    for i, di in enumerate(d):
        D=np.sum(np.sqrt(pi[:-di]*pi[di:]))
        N_mean=np.sum(np.abs(unbiased_rc_idx[:-di]-unbiased_rc_idx[di:])==1)
        N_mean/=len(unbiased_rc_idx)
        MU[i]=N_mean/D
        
    return MU
    
def sg_transmat(rc_bin, pi, MU, max_d=1):
    """
    Build the SGOOP transition rate matrix. 
    Attention: max_d>1 may result in less accuracy due to poor sampling of higher-order neighbor transitions.
    """
    d=np.arange(max_d)+1
    S=np.zeros((rc_bin, rc_bin))
    for i, di in enumerate(d):
        j=np.arange(rc_bin-di)
        k=j+di
        
        S[j,k] = MU[i]*np.sqrt(np.ma.divide(pi[j], pi[k]))
        S[k,j] = MU[i]*np.sqrt(np.ma.divide(pi[k], pi[j]))
    
    S[np.diag_indices(rc_bin)] = -np.sum(S, axis=0) # Diagonal terms
    
    return S
    
def eigenval(sg_transmat):
    """
    Compute eigenvalues and eigenvectors.
    """
    eigenValues, eigenVectors = np.linalg.eig(sg_transmat)
    idx = eigenValues.real.argsort()[::-1]     # Sorting by eigenvalues
    eigenValues = eigenValues[idx]        # Order eigenvalues
    eigenVectors = eigenVectors[:,idx]    # Order eigenvectors
    eigenExp = np.exp(eigenValues)        # Calculate exponentials
    
    return eigenValues.real, eigenExp.real, eigenVectors.real
    
def sgap(sg_transmat, nslow):
    """
    Calculate the spectral gap with self-consistent number of slow processes expected along the reaction coordinate.
    """
    E = eigenval(sg_transmat)
    
    return E[1][nslow-1]-E[1][nslow]
