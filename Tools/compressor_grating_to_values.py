import numpy as np

def compressor_grating_values(angle_of_incidence_deg=30.8, grating_line_density_mm=1280, grating_pair_displacement_mm=0, wavel_m=800e-9):
    """
    Calculates GDD and TOD values from the compressor grating position
    angle_of_incidence - degrees
    grating_line_density - l/mm
    grating_pair_separation - mm
    wavel - m

    returns GDD and TOD in s^2 and s^3
    """
    c = 299792458 # speed of light in m/s
    N = 2 # number of passes
    m = -1 # diffraction order
    d = 1/(grating_line_density_mm*1000) # grating period in m
    L = grating_pair_displacement_mm/1000 # separation in m
    theta = angle_of_incidence_deg*np.pi/180

    GDD = -((N*(m**2)*(wavel_m**3)*L)/(2*np.pi*(c**2)*(d**2))) * ((1 - ((-m*(wavel_m/d) - np.sin(theta))**2))**(-3/2))
    TOD = -((3*wavel_m)/(2*np.pi*c))*((1 + (wavel_m/d)*np.sin(theta) - (np.sin(theta)**2))/(1 - ((wavel_m/d) - np.sin(theta))**2)) * GDD
    return GDD, TOD

"""
GDD, TOD = compressor_grating_values(grating_pair_displacement_mm=3.937495614034377e-06*1000)
print(GDD*(10**30))
print(TOD*(10**45))

"""