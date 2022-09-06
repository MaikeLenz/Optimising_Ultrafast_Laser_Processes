import numpy as np
from scipy import integrate

def moment(x,y,n): 
    """
    returns an integral of the nth moment of intensity over frequency
    """
    integrand = np.array([]) #array of overlaps of the intensity with itself with a freq shift omega_shift
    for i in range(len(x)):
        integrand=np.append(integrand,y[i]*x[i]**n)
    return integrate.simps(integrand,x) 
   


def rms_width(x,y):
    """
    returns the rms width of an intensity distribution over angular frequencies omega
    """
    return 2 * np.sqrt(np.abs((moment(x,y,2)/moment(x,y,0))-(moment(x,y,1)/moment(x,y,0))**2))

"""
def pump_probe_width(omega, I):
    g = []

    for i in range(1, len(omega)+1):
        I_cut = I[:i]
        I_shifted = np.append(I_cut, np.zeros(len(omega)-i))
        g_tau = integrate.simps(I_shifted*I, omega)
        g.append(g_tau)
    for i in range(1, len(omega)):
        I_cut = I[i:]
        I_shifted = np.zeros(i)
        I_shifted = np.append(I_shifted, I_cut)
        g_tau = integrate.simps(I_shifted*I, omega)
        g.append(g_tau)
    
    tau_axis = np.append(-omega, omega[1:])
    
    N_2 = integrate.simps(g, tau_axis)
    variance_unnormalised = integrate.simps(g*(tau_axis**2), tau_axis)
    variance = variance_unnormalised/(2*N_2)
    print(variance)
    
def gauss(x, A, u, o):
    return A*np.exp(((x-u)**2)/(o**2))
"""