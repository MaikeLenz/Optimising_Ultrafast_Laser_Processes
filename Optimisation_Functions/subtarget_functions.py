# Target functions to be optimised for
##############################################################################################################################
# Imports
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np 
from scipy import integrate
import sys

sys.path.append('C:\\Users\\iammo\\Documents\\Optimising_Ultrafast_Laser_Processes\\Tools\\')
from rms_width import *
from get_phase import *
from envelopes import *
##############################################################################################################################
# Functions

def max_wavel_bandwidth(λ,Iλ):
    """
    Calculates the RMS bandwidth of a spectrum
    """
    return 10**9*rms_width(λ,Iλ)

def max_freq_bandwidth(λ,Iλ):
    """
    Converts a spectrum into the frequency domain and then calculates the RMS bandwidth
    """
    c = 299792458
    f = c/λ
    return rms_width(f,Iλ)

def min_duration_FT(om,Eom):
    """
    Fourier transforms to time domain and calculates the RMS width of the pulse
    """
    Et = np.fft.ifft(Eom)
    dom = om[2] - om[1]
    df = dom/(2*np.pi)
    t = np.fft.fftshift(np.fft.fftfreq(len(Et), d=df))
    return -rms_width(t,np.abs(Et)**2)

def min_thresh_duration_FT(om,Eom):
    """
    Fourier transforms to time domain to minimise duration
    """
    Et = np.fft.ifft(Eom)
    dom = om[2] - om[1]
    df = dom/(2*np.pi)
    t = np.fft.fftshift(np.fft.fftfreq(len(Et), d=df))
    return -threshold(t,np.abs(Et)**2)

def max_peak_power_FT(om,Eom):
    """
    Fourier transforms Eomega to Et and maximised amplitude there
    """
    Et = np.fft.ifft(Eom)
    #dom = om[2] - om[1]
    #df = dom/(2*np.pi)
    #t = np.fft.fftshift(np.fft.fftfreq(len(Et), d=df))

    #popt,_=curve_fit(gauss_envelope,t,np.abs(Et)**2, p0=[max(np.abs(Et)**2),2e-14,t[np.argmax(np.abs(Et)**2)]])
    #plt.plot(t,np.abs(Et)**2)
    #plt.plot(t,gauss_envelope(t,*popt))
    #plt.show()
    return max(np.abs(Et)**2)

def peak_power_window(λ,Iλ):
    """
    defined in the bossfunction
    """
    return True

def max_intens_integral(λ,Iλ,bounds):
    """
    integrates spectral intensity between wavelength bounds
    """
    λ_truncated1 = λ[bounds[0] <λ]
    λ_truncated=λ_truncated1[λ_truncated1< bounds[1]]
    Iλ_truncated1 = Iλ[bounds[0] <λ]
    Iλ_truncated=Iλ_truncated1[λ_truncated1< bounds[1]]
    return integrate.simps(Iλ_truncated,λ_truncated)

def threshold(x,y,x2=None,y2=None):
    """
    Find points where signal reaches certain level above noise and find the distance between them
    """
    thresh=0.1
    rows = np.where(y > max(y)*thresh)[0]
    if len(rows) <= 1:
        print("Set a lower threshold")
        min_index = rows[0]
        max_index = min_index
    else:
        min_index = rows[0]
        max_index = rows[-1]
    return x[max_index]-x[min_index]

def norm_and_int(x, y):
    """
    Normalises the pulse so that the maximum is 1, and then integrates
    """
    maximum = max(y)
    norm = []
    for i in y:
        norm.append(i/maximum)
    return integrate.simps(norm, x)

def max_peak_power_300nm(om,Eom):
    # First smooth using super Gaussian filter
    c=299792458
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 300e-9, 300e-9*0.1)
    #plt.plot(λ,np.abs(Eom)**2)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])
    #plt.plot(λ, np.abs(Eom_smooth)**2)
    #plt.show()
    # Now Fourier transform
    Et = np.fft.ifft(Eom_smooth)
    return max(np.abs(Et)**2)

def max_peak_power_300nm_envelope(om,Eom):
    # First smooth using super Gaussian filter
    c=299792458
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 300e-9, 300e-9*0.1)
    #plt.plot(λ,np.abs(Eom)**2)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])
    #plt.plot(λ, np.abs(Eom_smooth)**2)
    #plt.show()
    # Now Fourier transform
    Et = np.fft.ifft(Eom_smooth)
    dom = om[2] - om[1]
    df = dom/(2*np.pi)
    t = np.fft.fftshift(np.fft.fftfreq(len(Et), d=df))
    popt,_=curve_fit(gauss_envelope,t,np.abs(Et)**2, p0=[max(np.abs(Et)**2),2e-14,t[np.argmax(np.abs(Et)**2)]])

    return popt[0]

def max_peak_power_1300nm(om,Eom):
    # First smooth using super Gaussian filter
    c=299792458
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 1300e-9, 1300e-9*0.2)
    #plt.plot(λ,np.abs(Eom)**2)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])
    #plt.plot(λ, np.abs(Eom_smooth)**2)
    #plt.show()
    # Now Fourier transform
    Et = np.fft.ifft(Eom_smooth)
    return max(np.abs(Et)**2)

def max_peak_power_1200nm(om,Eom):
    # First smooth using super Gaussian filter
    c=299792458
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 1200e-9, 1200e-9*0.2)
    #plt.plot(λ,np.abs(Eom)**2)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])
    #plt.plot(λ, np.abs(Eom_smooth)**2)
    #plt.show()
    # Now Fourier transform
    Et = np.fft.ifft(Eom_smooth)
    return max(np.abs(Et)**2)
"""
def max_peak_power_300nm_quadratic_phase(om, Eom):
    # First get phase of pulse in freq domain
    om0 = moment(om,np.abs(Eom)**2,1)/moment(om,np.abs(Eom)**2,0) # Determine central frequency
    c = 299792458
    lambda0 = (2*np.pi*c)/om0
    phase = get_phase(om, Eom, lambda0)

    # Slice phase to only select part within pulse
    thresh = 0.1
    rows = np.where(np.abs(Eom)**2 > max(np.abs(Eom)**2)*thresh)[0]
    min_index = rows[0]
    max_index = rows[-1]
    phase_slice = phase[min_index-25:max_index+25]
    om_slice = om[min_index-25:max_index+25]

    # Fit a quadratic to the phase and determine the rms error
    def quad(x, a, b, c):
        return a*(x**2) + b*x + c
    quad_popt, _ = curve_fit(quad, om_slice, phase_slice, p0=[1,1,0])
    rms_phase_err = errorCorrection_int(om_slice, phase_slice, quad(om_slice, *quad_popt)) # Note: this value is negative

    # Smooth using super Gaussian filter
    c=299792458
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 300e-9, 300e-9*0.1)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])
    
    # Now Fourier transform
    Et = np.fft.ifft(Eom_smooth)
    dom = om[2] - om[1]
    df = dom/(2*np.pi)
    t = np.fft.fftshift(np.fft.fftfreq(len(Et), d=df))
    popt,_ = curve_fit(gauss_envelope,t,np.abs(Et)**2, p0=[max(np.abs(Et)**2),2e-14,t[np.argmax(np.abs(Et)**2)]])

    return popt[0] + 2*rms_phase_err
"""
def max_peak_power_300nm_quadratic_phase(om, Eom):
    # First get phase of pulse in freq domain
    om0 = moment(om,np.abs(Eom)**2,1)/moment(om,np.abs(Eom)**2,0) # Determine central frequency
    c=299792458
    lambda0 = (2*np.pi*c)/om0
    phase = get_phase(om, Eom, lambda0)

    # Smooth electric field using super Gaussian filter
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 300e-9, 300e-9*0.1)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])

    # Slice phase to only select part within pulse
    thresh = 0.1
    rows = np.where(np.abs(Eom_smooth)**2 > max(np.abs(Eom_smooth)**2)*thresh)[0]
    min_index = rows[0]
    max_index = rows[-1]
    phase_slice = phase[min_index-25:max_index+25]
    om_slice = om[min_index-25:max_index+25]

    # Fit a quadratic to the phase and remove this
    def quad(x, a, b, c):
        return a*(x**2) + b*x + c
    quad_popt, _ = curve_fit(quad, om_slice, phase_slice, p0=[1,1,0])
    phase_to_remove = quad(om_slice, *quad_popt)
    new_phase = np.zeros(len(om))
    for i in range(len(phase_to_remove)):
        new_phase[i+min_index-25] += phase_slice[i] - phase_to_remove[i]

    # Add the phase back to the intensity profile
    Eom_complex = []
    for i in range(len(om)):
        Eom_complex.append(np.abs(Eom_smooth[i])*np.exp(-1j*new_phase[i]))

    # Now Fourier transform
    Et = np.fft.ifft(Eom_complex)

    return max(np.abs(Et)**2)

def max_peak_power_1300nm_quadratic_phase(om, Eom):
    # First get phase of pulse in freq domain
    om0 = moment(om,np.abs(Eom)**2,1)/moment(om,np.abs(Eom)**2,0) # Determine central frequency
    c=299792458
    lambda0 = (2*np.pi*c)/om0
    phase = get_phase(om, Eom, lambda0)

    # Smooth electric field using super Gaussian filter
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 1300e-9, 1300e-9*0.2)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])

    # Slice phase to only select part within pulse
    thresh = 0.1
    rows = np.where(np.abs(Eom_smooth)**2 > max(np.abs(Eom_smooth)**2)*thresh)[0]

    try:
        min_index = rows[0]
        max_index = rows[-1]
        phase_slice = phase[min_index-25:max_index+25]
        om_slice = om[min_index-25:max_index+25]

        # Fit a quadratic to the phase and remove this
        def quad(x, a, b, c):
            return a*(x**2) + b*x + c
        quad_popt, _ = curve_fit(quad, om_slice, phase_slice, p0=[1,1,0])
        phase_to_remove = quad(om_slice, *quad_popt)
        new_phase = np.zeros(len(om))
        for i in range(len(phase_to_remove)):
            new_phase[i+min_index-25] += phase_slice[i] - phase_to_remove[i]

        # Add the phase back to the intensity profile
        Eom_complex = []
        for i in range(len(om)):
            Eom_complex.append(np.abs(Eom_smooth[i])*np.exp(-1j*new_phase[i]))

        # Now Fourier transform
        Et = np.fft.ifft(Eom_complex)

        return max(np.abs(Et)**2)
    except:
        return 0


def max_peak_power_1200nm_quadratic_phase(om, Eom):
    # First get phase of pulse in freq domain
    om0 = moment(om,np.abs(Eom)**2,1)/moment(om,np.abs(Eom)**2,0) # Determine central frequency
    c=299792458
    lambda0 = (2*np.pi*c)/om0
    phase = get_phase(om, Eom, lambda0)

    # Smooth electric field using super Gaussian filter
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 1200e-9, 1200e-9*0.2)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])
    
    try:
        # Slice phase to only select part within pulse
        thresh = 0.1
        rows = np.where(np.abs(Eom_smooth)**2 > max(np.abs(Eom_smooth)**2)*thresh)[0]
        min_index = rows[0]
        max_index = rows[-1]
        phase_slice = phase[min_index-25:max_index+25]
        om_slice = om[min_index-25:max_index+25]

        # Fit a quadratic to the phase and remove this
        def quad(x, a, b, c):
            return a*(x**2) + b*x + c
        quad_popt, _ = curve_fit(quad, om_slice, phase_slice, p0=[1,1,0])
        phase_to_remove = quad(om_slice, *quad_popt)
        new_phase = np.zeros(len(om))
        for i in range(len(phase_to_remove)):
            new_phase[i+min_index-25] += phase_slice[i] - phase_to_remove[i]

        # Add the phase back to the intensity profile
        Eom_complex = []
        for i in range(len(om)):
            Eom_complex.append(np.abs(Eom_smooth[i])*np.exp(-1j*new_phase[i]))

        # Now Fourier transform
        Et = np.fft.ifft(Eom_complex)

        return max(np.abs(Et)**2)
    except:
        return 0

def max_peak_power_200nm_quadratic_phase(om, Eom):
    # First get phase of pulse in freq domain
    om0 = moment(om,np.abs(Eom)**2,1)/moment(om,np.abs(Eom)**2,0) # Determine central frequency
    c=299792458
    lambda0 = (2*np.pi*c)/om0
    phase = get_phase(om, Eom, lambda0)

    # Smooth electric field using super Gaussian filter
    λ=(2*np.pi*c)/om
    filter = superGauss(λ, 200e-9, 200e-9*0.1)
    Eom_smooth = []
    for i in range(len(Eom)):
        Eom_smooth.append(Eom[i]*filter[i])
    
    try:
        # Slice phase to only select part within pulse
        thresh = 0.1
        rows = np.where(np.abs(Eom_smooth)**2 > max(np.abs(Eom_smooth)**2)*thresh)[0]
        min_index = rows[0]
        max_index = rows[-1]
        phase_slice = phase[min_index-25:max_index+25]
        om_slice = om[min_index-25:max_index+25]

        # Fit a quadratic to the phase and remove this
        def quad(x, a, b, c):
            return a*(x**2) + b*x + c
        quad_popt, _ = curve_fit(quad, om_slice, phase_slice, p0=[1,1,0])
        phase_to_remove = quad(om_slice, *quad_popt)
        new_phase = np.zeros(len(om))
        for i in range(len(phase_to_remove)):
            new_phase[i+min_index-25] += phase_slice[i] - phase_to_remove[i]

        # Add the phase back to the intensity profile
        Eom_complex = []
        for i in range(len(om)):
            Eom_complex.append(np.abs(Eom_smooth[i])*np.exp(-1j*new_phase[i]))

        # Now Fourier transform
        Et = np.fft.ifft(Eom_complex)

        return max(np.abs(Et)**2)
    except:
        return 0