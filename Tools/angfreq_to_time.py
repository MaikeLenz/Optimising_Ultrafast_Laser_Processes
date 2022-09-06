#takes complex e(omega) as input and returns E(t)

#from locale import ERA_T_FMT
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.fft import *

def E_omega_to_E_time(E_om,om):
    """
    Fourier transforms from ang freq to time domain 
    """
    Fs=2*(om[np.nonzero(E_om)[-1]]) #samplicng frequency
    E_t = sp.fft.ifft(E_om)
    print(E_t)
    t = np.linspace(0,len(E_t)/Fs,len(E_t))
    #t=np.append(t[len(t)/2:],t[:len(t)/2])
    E_t=np.append(E_t[int(len(E_t)/2):],E_t[:int(len(E_t)/2)])
    plt.plot(t,E_t)
    plt.show()
    return E_t,t

def FourierTransform(Et,t):
    E=sp.fft.fft(Et)
    omf=sp.fft.fftfreq(len(t))
    """
    #samp=(abs(abs(x[-1])-abs(x[0]))/len(x))/1000 
    tsamp=(abs(abs(t[-1])-abs(t[0]))/len(t))/1000
    om=omf[int(len(omf)/2+1):len(omf)]
    #repx= 2*tsamp/om
    y = abs(E[int(len(omf)/2+1):len(omf)]).tolist()
    """
    plt.plot(omf,E)
    plt.show()

def t_to_f(t, Et):
    """
    Input of t and Et
    Returns f, Ef
    """
    """
    ts = t[1]-t[0] # Sampling interval
    sr = 1/ts # Sampling rate
    Ef = fft(Et)
    N = len(Ef) # Number of points
    n = np.arange(N)
    T = N/sr
    f = n/T

    Ef_oneside = Ef[:N//2]
    f_oneside = f[:N//2]
    return f_oneside, Ef_oneside
    """
    Ef = fft(Et)
    Ef = fftshift(Ef) # shift zero-frequency component to center of spectrum

    N = len(Et) # Number of points

   
    dt=t[1]-t[0]
    #F=1/dt
    #df=F/len(Et)
    #f=np.arange(0,F,df)
    f = fftfreq(N, dt)
    f = fftshift(f) # shift zero-frequency component to center of spectrum

    #Ef_oneside = list(Ef[:N//2])
    #Ef_otherside=list(Ef[N//2:])
    #Ef_new=np.array(Ef_otherside+Ef_oneside)
    #t_oneside = t[:N//2]
    return f,Ef

def f_to_t(f, Ef):
    """
    Input of f and Ef
    Returns t, Et
    """
    if any(i<0 for i in f) == False:
        Ef=np.append(Ef, Ef)
    else:
        neg_indices=np.where(f < 0)[0]
        pos_indices=np.where(f>=0)[0]
        Ef_neg=Ef[neg_indices]
        Ef_pos=Ef[pos_indices]
        Ef=np.array(list(Ef_pos)+list(Ef_neg))
        #need to put these values at the end of the Et array s.t. first half is positive and sceond is -ve frequencies
    Et = ifft(Ef)
    N = len(Ef) # Number of points
    #t = 1/f
    """
    fs = np.abs(f[1]-f[0])
    sr = 1/fs
    N = len(Et)
    n = np.arange(N)
    T = N/sr
    t = n/T
    """
    """
    #sampling frequency is Fs
    Fs=2*max(f)
    t = np.arange(0, (N-1)/Fs, 1/Fs)
    """

def f_to_t_irfft(f, Ef):
    """
    Input of f and Ef
    Returns t, Et
    """
    Et = np.fft.irfft(Ef)
    N = len(Ef) # Number of points
    #t = 1/f
    df=np.abs(f[1]-f[0])#smallest frequency difference gives inverse of duration
    T=1/df#overall duration
    print("T=",T)
    t=np.linspace(-T/2,T/2,len(Et))#construct time axis
    print("t:",t[0],t[-1])
    Et_oneside = list(Et[:N//2])
    Et_otherside=list(Et[N//2:])
    Et_new=np.array(Et_otherside+Et_oneside)

    return t,Et
"""
t = np.linspace(0,1,100)
freq = 2
Et = 3*np.sin(2*np.pi*freq*t)
plt.plot(t, Et, label='Before FFT')
plt.legend()

f, Ef = t_to_f(t, Et)
plt.figure()
plt.plot(f, Ef, label='After FFT')
plt.legend()

t2, Et2 = f_to_t(f, Ef)
plt.figure()
plt.plot(t2, Et2, label='After IFFT')
plt.legend()
plt.show()
"""
"""
from scipy.fft import fft, fftfreq, fftshift
def get_intensity_spectrum(t, Et):

#    Uses scipy.fft to calculate the intensity spectrum of a laser pulse E(t). 
#    `t` is time axis, `Et` is (complex) laser electric field sampled on t
#    returns tuple (`omega`, `I`), where `omega` is angular frequency and `I` is spectral intensity.
    
#    Tip: use 1024 or 2048 time points

    assert len(t) == len(Et)
 
    t = np.array(t)
    Et = np.array(Et)
    
    N = len(t) # num points
    dt = t[1]-t[0] # time step
    f = fftfreq(N, dt) # frequency axis
    f = fftshift(f) # shift zero-frequency component to center of spectrum
    omega = 2 * np.pi * f # angular frequency axis

    Ef = fft(Et) # the fft
    Ef = fftshift(Ef) # shift zero-frequency component to center of spectrum
    I = np.abs(Ef)**2
    I /= max(I)
   
    return omega, I
"""