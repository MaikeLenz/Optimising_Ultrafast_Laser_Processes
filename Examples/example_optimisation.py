import sys
import csv

path_to_repo = 'C:\\Users\\iammo\\Documents\\Optimising_Ultrafast_Laser_Processes\\'
#path_to_repo = 'D:\\HiDrive\\users\\maikelenz\\docs\\MSci_FinalCode\\Optimising_Ultrafast_Laser_Processes\\'

sys.path.append(path_to_repo+'Optimisation_Functions\\')
from optimisation_function import *

# First choose the paramters you wish to vary
# Values:  radius, flength, gas, pressure, wavelength, energy, τfwhm, grating_pair_separation
params=["energy", "pressure", "grating_pair_displacement", "radius", "flength"]

# Define experimental parameters
# Initial values must be chosen even for the parameters you are optimising so that an initial simulation can be run
gas="He" # The gas to fill the HCF with
radius = 100e-6 # HCF core radius
flength = 3 # HCF length
wavel = 800e-9 # Central wavelength of the input laser pulse
duration = 30e-15 # Initial pulse duration
energy_init = 0.75e-3 # Energy of the input pulse
pressure_init = 0.66*3 # Pressure within the HCF
grating_init = 0 # Compressor grating separation

initial_values_HCF = [radius, flength, gas, pressure_init, wavel, energy_init, duration, grating_init] # These values must remain in the order given here to be passed into the optimisation function

# You can also choose the bounds of the parameters to vary and specify these in a dictionary
# If no bounds are chosen, the optimisation will use the default parameter bounds
bounds = {"energy":(0.5e-3, 1.5e-3), "pressure":(0.66*1.0,8.0*0.66), "grating_pair_displacement":(-0.5e-3, 0.5e-3), "radius":(50e-6, 500e-6), "flength":(0.1, 10)}

# Choose the number of initial points and number of iterations for Bayesian optimisation to run for
# We would recommend starting with one initial point and one iteration to test the optimisation function
inits = 1
iters = 1

# Run the optimisation function
result, iterations = Luna_BO(params, initial_values_HCF, function=max_peak_power_1300nm_quadratic_phase, 
                            init_points=inits, n_iter=iters, parameter_bounds=bounds, subtarget_analysis="f")
# result contains only the final optimal parameters. iterations contains the optimal parameters found after every iteration of the optimisation function.

# Extract results
target = result['target']
grating = result["params"]["grating_pair_displacement"]
pressure = result["params"]["pressure"]
energy = result["params"]["energy"]
radius = result["params"]["radius"]
flength = result["params"]["flength"]

# Below is a suggestion for how the final optimal parameters can be saved in a csv file
savefile_name = path_to_repo+'example_optimisation_results.csv'
header = ['init_points', 'n_iter', 'peak power', 'energy, J', 'pressure, bar', 'radius, m', 'flength, m', 'FWHM, s', 'wavel, m', 'gas', 'grating_pair_displacement, m']
with open(savefile_name, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    # write the header
    writer.writerow(header)
    # write the data
    writer.writerow([inits, iters, target, energy, pressure, radius, flength, duration, wavel, gas, grating])
