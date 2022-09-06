import sys
import csv
#sys.path.append('C:\\Users\\ML\\OneDrive - Imperial College London\\MSci_Project\\code\\Synth\\Optimising-Field-Synthesiser\\HCF sim\\Python\\Luna_BO\\')
sys.path.append('C:\\Users\\iammo\\Documents\\Optimising-Field-Synthesiser\\HCF sim\\Python\\Luna_BO\\')
from bossfunction_Luna_debugging_Mollypaths import *

params=["energy", "pressure", "grating_pair_displacement", "radius", "flength"]

#values:  radius, flength, gas, pressure, wavelength, energy, τfwhm, grating_pair_separation
gas="He"
# Define experimental params

radius = 100e-6 # HCF core radius
flength = 3 # HCF length
wavel=800e-9
duration=30e-15
energy_init=0.75e-3
pressure_init=0.66*3
grating_init=0

initial_values_HCF=[radius, flength, gas,pressure_init, wavel, energy_init,duration, grating_init]


inits=1
iters=1
result,iterations=Luna_BO_debug(params, initial_values_HCF, function=max_peak_power_1300nm_quadratic_phase, init_points=inits, n_iter=iters)
target = result['target']
grating=result["params"]["grating_pair_displacement"]
pressure=result["params"]["pressure"]
energy=result["params"]["energy"]
radius=result["params"]["radius"]
flength=result["params"]["flength"]

"""
header = ['init_points', 'n_iter', 'peak power', 'energy, J', 'pressure, bar', 'radius, m', 'flength, m', 'FWHM, s', 'wavel, m', 'gas', 'grating_pair_displacement, m']
#with open('C:\\Users\\iammo\\Documents\\Optimising-Field-Synthesiser\\HCF sim\\Python\\Luna_BO\\tests\\optimise_Luna\\data\\peak_power_1200e-9wavelwindow_varyfwhm__init_' + str(init_points) + '_niter_' + str(n_iter) + '.csv', 'a', encoding='UTF8', newline='') as f:
with open('C:\\Users\\iammo\\Documents\\Optimising-Field-Synthesiser\\HCF sim\\Python\\Luna_BO\\tests\\optimise_with_FT_and_phase\\data\\1300nm_phasecondition\\1300nm_He__init_50_niter_1000.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f) # peak_power_1000e-9wavelwindow__init_50_niter_100
    # write the header
    writer.writerow(header)

    # write the data
    writer.writerow([inits, iters, target, energy, pressure, radius, flength, duration, wavel, gas, grating])
"""
