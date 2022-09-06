# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 11:55:57 2019

@author: jt
"""

import numpy as np
import matplotlib.pylab as plt

def efield_time_domain(t, amp, om0, dom, t0, gdd, cep):
    """ returns time domain e_efield (real)
    t is time axis
    amp is pulse amplitude
    om0 is centre frequency
    dom is bandwidth
    t0 is time offset
    gdd is group delay dispersion (in units of time^2)
    cep is carrier envelop phase
    """
    
    q = (np.log(4)/dom **2 + 0.5 * 1j * gdd )**0.5
    return np.real(amp * np.log(4)**0.5 * np.exp(-0.25*(t-t0)**2/q**2) * np.exp(1j*om0*(t-t0) + 1j*cep)/(dom*q))


def efield_freq_domain(t, amp, om0, dom, t0, gdd, cep):
    """ returns freq domain e_field (real) by fourier transforming time domain
    t is time axis
    amp is pulse amplitude
    om0 is centre frequency
    dom is bandwidth
    t0 is time offset
    gdd is group delay dispersion (in units of time^2)
    cep is carrier envelop phase
    """
    
    q = ( np.log(4)/dom **2 + 0.5 * 1j * gdd )**0.5
    E=amp * np.log(4)**0.5 * np.exp(-0.25*(t-t0)**2/q**2) * np.exp(1j*om0*(t-t0) + 1j*cep)/(dom*q)
    E_omega=np.fft.ifft(E)
    omega=np.fft.fftfreq(len(t),d=(t[1]-t[0]))
    return np.real(E_omega), omega

def get_phi(omega, omega0, CEP, GD, GDD, TOD, FoOD, FiOD):
    return CEP + GD*(omega-omega0) + (1/2)*GDD*(omega-omega0)**2 + (1/6)*TOD*(omega-omega0)**3 + (1/24)*FoOD*(omega-omega0)**4 + (1/120)*FiOD*(omega-omega0)**5

def E_field_freq(omega, GD=0.0, wavel=1000, domega=2, amp=1, CEP=0, GDD=0, TOD=0, FoOD=0, FiOD=0):
    """
    Defines an E-field pulse in the frequency domain
    omega: array of angular frequencies
    GD: group delay (relates to absolute arrival time t0)
    wavel: wavelength
    domega: bandwidth
    amp: amplitude
    CEP: carrier envelope phase
    GDD: group delay dispersion
    TOD: third order dispersion
    """
    c = 299792458 # m/s
    E0 = amp
    omega0 = 2 * np. pi * c/wavel # rad/fs
    E_transform_limited = E0 * np.exp(-2 * np.log(2) * (omega-omega0)**2/domega**2)
    phi = get_phi(omega, omega0, CEP, GD, GDD, TOD, FoOD, FiOD)

    E = E_transform_limited * np.exp(phi * 1J)
    phase_wrapped = np.angle(E)
    phase = np.unwrap(phase_wrapped)
    return E, phase