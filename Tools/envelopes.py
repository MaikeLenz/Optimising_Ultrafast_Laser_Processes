import numpy as np

def gauss_envelope(t,A,sig,mu):
    """
    Gaussian envelope
    """
    return A*np.exp(-(t-mu)**2/sig**2)

def superGauss(x,x0,sigma):
    """
    Super Gaussian envelope
    """
    return np.exp(-((x-x0)/sigma)**4)
