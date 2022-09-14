# The Optimisation Function

Optimisations are done by calling the "Luna_BO" function in the optimisation_function.py file. 

This takes the following arguments:
- params: A list of strings. This denotes the parameters that are going to be optimised. Options are 
  - "radius". Core radius of the HCF in m.
  - "flength". Fibre length in m.
  - gas type. e.g. "He", "Ne", "Ar"
  - "pressure". The average gas pressure in the fibre in Bar.
  - "λ0". The central wavelength of the input laser pulse in m.
  - "energy". Energy of the input laser pulse in J.
  - "FWHM". FWHM duration of the input laser pulse in s.
  - "grating_pair_displacement". Compressor grating plate separation in m.
