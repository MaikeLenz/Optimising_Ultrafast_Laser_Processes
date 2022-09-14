# The Optimisation Function

Optimisations are done by calling the "Luna_BO" function in the optimisation_function.py file. 

This takes the following arguments:
- **params**: A list of strings. This denotes the parameters that are going to be optimised. Options are 
  - "radius". Core radius of the HCF in m.
  - "flength". Fibre length in m.
  - gas type. e.g. "He", "Ne", "Ar"
  - "pressure". The average gas pressure in the fibre in bar.
  - "λ0". The central wavelength of the input laser pulse in m.
  - "energy". Energy of the input laser pulse in J.
  - "FWHM". FWHM duration of the input laser pulse in s.
  - "grating_pair_displacement". Compressor grating plate separation in m.
- **initial_values_HCF**: The starting point of the Optimisation. This is a list of the following values: [radius, flength, gas_type, pressure_init, wavel, energy, FWHM, grating_pair_displacement] It must be given in this particular order!
- **function**. This is the function that is carried out in every iteration to return the value to be optimised. Select one from the subtarget_functions.py file
- **init_points**. The number of initial random searches in parameter that are carried out before optimising. This is optional; the default is set to 50.
- **n_iter**. The number of iterations of the optimisation that are going to be carried out. This is optional; teh default is set to 50.
- **subtarget_analysis**. This determines whether the subtarget function will be carried out in terms of wavelength or frequency. It can be set to either "f" for frequency or "w" for wavelength. This is optional; the default is frequency.
- **Gaussian**. This is a boolean. If false, a custom input electric field is defined in frequency with spectral phase. This is found from the compressor grating separation. If true, a default Gaussian pulse shape in frequency is used.
- **ImperialLab**. This is a boolean. If true, experimentally measured input spectra (in the file Imperial_Lab_Data\Input_Power_Scan.txt) with spectral phase as calculated from the compressor grating separation are passed into Luna. If false, the spectral spread is calculated from the FWHM duration of the pulse and spectral phase is again calculated from the compressor grating separation.
- **parameter_bounds**. This is a boolean. If false, default parameter bounds are used which are defined in optimisation_function.py. If true, the selected bounds are used instead.

Inside the Luna_BO function, a native target function is defined which will return the value to be optimised. This is essentially just the subtarget function that was selected but with the correct arguments.


Luna_BO returns:
- **Results**. This is a dictionary of the final parameters which can be accessed as e.g. results["params"]["grating_pair_displacement"] or
results["params"]["pressure"] and also the final output of the subtarget function which can be accessed as results['target'].
- **Iterations**. This every point in parameter space that was probed by the optimisation.


