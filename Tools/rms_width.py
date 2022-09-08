import numpy as np
from scipy import integrate

def moment(x,y,n): 
    """
    Returns an integral of the nth moment of intensity over frequency
    """
    integrand = np.array([]) #array of overlaps of the intensity with itself with a freq shift omega_shift
    for i in range(len(x)):
        integrand=np.append(integrand,y[i]*x[i]**n)
    return integrate.simps(integrand,x) 

def rms_width(x,y):
    """
    Returns the RMS width of an intensity distribution over angular frequencies omega
    """
    return 2 * np.sqrt(np.abs((moment(x,y,2)/moment(x,y,0))-(moment(x,y,1)/moment(x,y,0))**2))
