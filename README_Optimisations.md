# The Optimisation Function

Optimisations are done by calling the "Luna_BO" function in the optimisation_function.py file. 

This takes the following arguments:
    * params: A list of strings. This denotes the parameters that are going to be optimised. Options are 
        * "radius". Core radius of the HCF
        * "flength". Fibre length
        * the gas type. e.g. "He", "Ne", "Ar"
        * "pressure". The average gas pressure in the fibre. 
        * "λ0". The central wavelength of the input laser pulse 
        * "energy". Energy of the input laser pulse 
        * "FWHM". FWHM duration of the input laser pulse 
        * "grating_pair_displacement"