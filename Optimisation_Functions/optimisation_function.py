# This function carries out Bayesian Optimisation for hollow core fibres
#params to be varied: 
    # Pulse: input energy, τfwhm, central wavelength
    # Fibre: pressure, fibre core radius, fibre length
##############################################################################################################################
# Imports
import numpy as np
from bayes_opt import BayesianOptimization
import matplotlib.pyplot as plt
import csv
import julia

# Add your Julia path here
julia_path="C:\\Users\\iammo\\AppData\\Local\\Programs\\Julia-1.7.1\\bin\\julia.exe"
#julia_path="C:\\Users\\ML\\AppData\\Local\\Programs\\Julia-1.7.0\\bin\\julia.exe"

julia.Julia(runtime=julia_path)
from julia import Main

import sys

#path to github repository
path_to_repo='C:\\Users\\iammo\\Documents\\Optimising_Ultrafast_Laser_Processes\\'
#path_to_repo='D:\\HiDrive\\users\\maikelenz\\docs\\MSci_FinalCode\\Optimising_Ultrafast_Laser_Processes\\'

#import necessary functions
sys.path.append(path_to_repo+'Optimisation_Functions\\')
sys.path.append(path_to_repo+'Tools\\')
from pulse_with_GDD import *
from subtarget_functions import *
from compressor_grating_to_values import *

##############################################################################################################################
c = 299792458 # m/s

#here starts th eoptimisation function
def Luna_BO(params, initial_values_HCF, function, init_points = 50, n_iter = 50, subtarget_analysis = "f", 
    Gaussian = False, ImperialLab = False, parameter_bounds = None, 
    plotting = True, wavel_bounds = None,
    record_iterations = False, save_path = None):     
    """
    performs BO with params as specified as strings in params input (params is list of strings) on the HCF.
    init_points: number of initial BO points
    n_iter: number of iterations
    plots input&output spectrum
    initial_values_HCF = [radius, flength, gas, pressure, λ0, τfwhm, energy] # array of initial values for Luna simulation, could change this to input actual Luna simulation
    """ 
    # Start by assigning values to Luna simulation
    Main.using("Luna")
    
    Main.radius = initial_values_HCF[0]
    Main.flength = initial_values_HCF[1]
    Main.gas_str = initial_values_HCF[2]
    Main.eval("gas = Symbol(gas_str)")
    Main.pressure = initial_values_HCF[3]
    Main.λ0 = initial_values_HCF[4]
    Main.energy = initial_values_HCF[5]
    Main.τfwhm = initial_values_HCF[6]
    grating_pair_displacement = initial_values_HCF[7]

    args_BO = {} #this dictionary will contain only the parameters we want to vary here
    params_dict={}
    params_dict['radius'] = initial_values_HCF[0]
    params_dict['flength'] = initial_values_HCF[1]
    params_dict['gas_str'] = initial_values_HCF[2]
    params_dict['pressure'] = initial_values_HCF[3]
    params_dict['λ0'] = initial_values_HCF[4]
    params_dict['energy'] = initial_values_HCF[5]
    params_dict['FWHM'] = initial_values_HCF[6]
    params_dict['grating_pair_displacement'] = initial_values_HCF[7]

    for i in params:
        if i in params_dict:
            args_BO[i] = params_dict[i] #append parameters to be varied to dictionary
    print(args_BO)

    def target_func(**args):
        """
        this is the target function of the optimiser. It is created as a nested function to take only the desired variables as inputs.
        It will consist of one of the sub-target functions in the subtarget function file or one of the rms error functions in ErrorCorrection_integrate.
        """
        for i in range(len(params)):
            args_BO[params[i]] = args[params[i]]
            params_dict[params[i]]=args[params[i]]
            
        # Update the simulation's variables with new parameters
        for key, value in args_BO.items():
            if 'energy' in key:
                Main.energy = value
            elif 'λ0' in key:
                Main.λ0 = value
            elif 'pressure' in key:
                Main.pressure = value
            elif 'radius' in key:
                Main.radius = value
            elif 'flength' in key:
                Main.flength = value
            elif 'FWHM' in key:
                Main.τfwhm = value
            elif 'grating_pair_displacement' in key:
                grating_pair_displacement = value

        # Critical power condition
        Main.eval('ω = PhysData.wlfreq(λ0)')
        Main.eval('_, n0, n2  = Tools.getN0n0n2(ω, gas; P=pressure)')
        Main.eval('Pcrit = Tools.Pcr(ω, n0, n2)')
        Pcrit = Main.Pcrit
        Pmin = 0
        τfwhm = Main.τfwhm
        tau = τfwhm/(2*np.sqrt(np.log(2)))
        P = Main.energy/(np.sqrt(np.pi)*tau)
        power_condition = int(Pmin <= P <= Pcrit)

    
        if Gaussian == False:
            if ImperialLab == False:
                # Custom data pulse is defined and passed to prop capillary
                domega = 2*np.pi*0.44/τfwhm
                c=299792458
                omega = np.linspace(2*np.pi*c/params_dict["λ0"] - 5*domega/2, 2*np.pi*c/params_dict["λ0"] + 5*domega/2, 1000)

                GDD, TOD = compressor_grating_values(grating_pair_displacement_mm=params_dict["grating_pair_displacement"]*1000)

                E, ϕω = E_field_freq(omega, GD=0.0, wavel=params_dict["λ0"], domega=domega, amp=1, CEP=0, GDD=GDD, TOD=TOD)
                Iω = np.abs(E)**2


                Main.ω = omega
                Main.Iω = Iω  
                Main.phase = ϕω 
                # Pass data to Luna
                Main.eval('pulse = Pulses.DataPulse(ω, Iω, phase; energy, λ0=NaN, mode=:lowest, polarisation=:linear, propagator=nothing)')
                Main.duv = Main.eval('duv = prop_capillary(radius, flength, gas, pressure; λ0, pulses=pulse, trange=400e-15, λlims=(150e-9, 4e-6))')
            else:
                # Get experimental input spectrum, scale energy and apply phase
                lines=[]
                columns=[[],[],[],[],[],[],[],[],[],[],[]]
                ### NEED TO FIX PATH IN LINE BELOW
                with open (path_to_repo+'Imperial_Lab_Data\\Input_Power_Scan.txt', 'rt') as myfile:  # Open lorem.txt for reading
                    for myline in myfile:
                        lines.append(myline)

                c = 299792458 # m/s
                data=lines[22:] #gets rid of all the stuff at the top
                data=data[int(len(data)/2):]
                for i in data:
                    cut=i.split("\t") #delimiter is \t
                    for j ,value in enumerate(cut):
                        columns[j].append(float(value))

                wavel_nm = np.array(columns[0])
                intens = np.array(columns[3])
                omega_list=2*np.pi*c/(wavel_nm*10**-9)
                Main.ω = omega_list[::-1]
                Main.Iω = intens[::-1]

                GDD, TOD = compressor_grating_values(grating_pair_displacement_mm=params_dict["grating_pair_displacement"]*1000)
                phase = []
                for j in range(len(omega_list)):
                    phase.append(get_phi(omega=omega_list[j], omega0=2*np.pi*c/params_dict["λ0"], CEP=0, GD=0, GDD=GDD, TOD=TOD, FoOD=0, FiOD=0))
                Main.phase = phase

                Main.eval('pulse = Pulses.DataPulse(ω, Iω, phase; energy, λ0=NaN, mode=:lowest, polarisation=:linear, propagator=nothing)')
                Main.duv = Main.eval('duv = prop_capillary(radius, flength, gas, pressure; λ0, pulses=pulse, trange=400e-15, λlims=(150e-9, 4e-6))')
                
        else:
            # Default Gaussian pulse passed to prop capillary
            Main.duv = Main.eval('duv = prop_capillary(radius, flength, gas, pressure; λ0, τfwhm, energy, trange=400e-15, λlims=(150e-9, 4e-6))')

        Main.eval("λ, Iλ = Processing.getIω(duv, :λ, flength)")
        Main.eval("grid = Processing.makegrid(duv)")
        Main.eval("ω = grid.ω")
        Main.eval('Eω = duv["Eω"][:,end]')
        """
        if function == peak_power_window:
            print("Evaluating peak power in wavelength bounds")
            Main.wavel_bounds=wavel_bounds
            #can also have Processing.peakpower
            Main.eval("peak_power=Stats.peakpower(duv, λlims=wavel_bounds)")
            peak_power=Main.peak_power
            print(peak_power)
            return peak_power*power_condition
        """
        # Get values    
        λ = Main.λ
        Iλ = Main.Iλ
        Iλ=Iλ.reshape(len(Iλ),)
        omega = Main.ω
        Eomega = Main.Eω  

        # Pass output from Luna into the subtarget function
        if subtarget_analysis == "w":
            return function(λ, Iλ, wavel_bounds)*power_condition
        elif subtarget_analysis == "f":
            return function(omega, Eomega)*power_condition
        else:
            raise ValueError("subtarget_analysis value is invalid. Must choose 'f' for frequency analysis or 'w' for wavelength analysis.")

    # Make pbounds dictionary 
    if parameter_bounds == None:
        pbounds = {}
        for i in params:
            #assume standard bounds
            if 'energy' in i:
                pbounds[i] = (0.5e-3, 1.5e-3)
            elif 'FWHM' in i:
                pbounds[i] = (20e-15, 35e-15)
            elif 'λ0' in i:
                pbounds[i] = (700e-9,900e-9)
            elif 'pressure' in i:
                if params_dict["gas_str"]=="He":
                    pbounds[i]=(0.66*1.0,8.0*0.66)
                if params_dict['gas_str']=="Ar":
                    pbounds[i] = (0.66*0.6, 0.66*1.5)
                elif params_dict['gas_str']=="Ne":
                    pbounds[i] = (0.66*3.0, 0.66*3.5)
                else:
                    pbounds[i]=(0.66*1.0,8.0*0.66)
            elif 'radius' in i:                
                pbounds[i] = (50e-6, 500e-6)
            elif 'flength' in i:
                pbounds[i] = (0.1, 10)
            elif 'grating_pair_displacement' in i:
                pbounds[i] = (-0.5e-3, 0.5e-3)
    else:
        pbounds = parameter_bounds

    print(pbounds)

    # Now create BO with the defined target function
    optimizer = BayesianOptimization(
        f=target_func,
        pbounds=pbounds,
        verbose=2, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state=1,
        )

    # Probe theoretically predicted SPM optimum
    if params_dict['gas_str']=="Ar":
        optimizer.probe(params={"energy": 1.1e-3, "pressure": 0.66*1.0, "grating_pair_displacement":0.0},lazy=True,)
        #optimizer.probe(params={"grating_pair_displacement":0.0},lazy=True,)
    elif params_dict['gas_str']=="Ne":
        optimizer.probe(params={"energy": 1.1e-3, "pressure": 0.66*3.5, "grating_pair_displacement":0.0},lazy=True,)
        #optimizer.probe(params={"grating_pair_displacement":0.0},lazy=True,)

    if record_iterations == True:
        # If True, csv files will be saved with the optimums during the optimisation, and the points probed during the initial random search and the optimisation
        # First probe the initial points
        optimizer.maximize(
            init_points=init_points,
            n_iter=0,
            )
        columns = ['iteration', 'target']
        for key, value in args_BO.items():
            columns.append(key)
        with open(save_path + '_optimums.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(columns)
        
        # Now probe each iteration and record the results
        for i in range(n_iter):
            optimizer.maximize(
            init_points=0,
            n_iter=1,
            )
            result = optimizer.max
            with open(save_path + '_optimums.csv', 'a', encoding='UTF8', newline='') as f:
                writer = csv.writer(f)
                row = [i, result["target"]]
                for key, value in args_BO.items():
                    row.append(result["params"][key])
                writer.writerow(row)

        # Record all points probed
        iterations = optimizer.res
        with open(save_path + '_points_probed.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            for i, res in enumerate(iterations):
                row = [i, res["target"]]
                for key, value in args_BO.items():
                    row.append(res["params"][key])
                writer.writerow(row)
    else:
        optimizer.maximize(
            #maximises the target function output. In the case of the RMS error functions, this is a minimisation because the errors are multiuplied by -1
            init_points=init_points,
            n_iter=n_iter,
            #acq="ucb", 
            #kappa=0.1
            )

    print(optimizer.max) #final parameters
    results=optimizer.max["params"]
    for key, value in args_BO.items():
        if 'energy' in key:
            Main.energy = results[key]
        elif 'λ0' in key:
            Main.λ0 = results[key]
        elif 'pressure' in key:
            Main.pressure = results[key]
        elif 'radius' in key:
            Main.radius = results[key]
        elif 'flength' in key:
            Main.flength = results[key]
        elif 'FWHM' in key:
            Main.τfwhm = results[key]   
        elif 'grating_pair_separation' in key:
            grating_pair_displacement = results[key]

    if plotting == True:
        if Gaussian == False:
            if ImperialLab == False:
                λ0 = Main.λ0
                τfwhm = Main.τfwhm
                domega = 2*np.pi*0.44/τfwhm
                c=299792458#m/s
                omega = np.linspace(2*np.pi*c/λ0 - 5*domega/2, 2*np.pi*c/λ0 + 5*domega/2, 100)

                GDD, TOD = compressor_grating_values(grating_pair_displacement_mm=grating_pair_displacement*1000)

                E, ϕω = E_field_freq(omega, GD=0.0, wavel=λ0, domega=domega, amp=1, CEP=0, GDD=GDD, TOD=TOD)
                Iω = np.abs(E)**2

                Main.ω = omega
                Main.Iω = Iω  
                Main.phase = ϕω 
                Main.eval('pulse = Pulses.DataPulse(ω, Iω, phase; energy, λ0=NaN, mode=:lowest, polarisation=:linear, propagator=nothing)')
                Main.duv = Main.eval('duv = prop_capillary(radius, flength, gas, pressure; λ0, pulses=pulse, trange=400e-15, λlims=(150e-9, 4e-6))')
            else:
                # Get experimental input spectrum, scale energy and apply phase
                lines=[]
                columns=[[],[],[],[],[],[],[],[],[],[],[]]
                with open (path_to_repo+'Imperial_Lab_Data\\Input_Power_Scan.txt', 'rt') as myfile:  # Open lorem.txt for reading
                    for myline in myfile:
                        lines.append(myline)

                c = 299792458 # m/s
                #print(lines[:22])
                data=lines[22:] #gets rid of all the stuff at the top
                data=data[int(len(data)/2):]
                for i in data:
                    cut=i.split("\t") #delimiter is \t
                    for j ,value in enumerate(cut):
                        columns[j].append(float(value))

                wavel_nm = np.array(columns[0])
                intens = np.array(columns[3])
                omega_list=2*np.pi*c/(wavel_nm*10**-9)
                Main.ω = omega_list[::-1]
                Main.Iω = intens[::-1]

                GDD, TOD = compressor_grating_values(grating_pair_displacement_mm=grating_pair_displacement*1000)
                phase = []
                for j in range(len(omega_list)):
                    phase.append(get_phi(omega=omega_list[j], omega0=2*np.pi*c/params_dict["λ0"], CEP=0, GD=0, GDD=GDD, TOD=TOD, FoOD=0, FiOD=0))
                Main.phase = phase

                Main.eval('pulse = Pulses.DataPulse(ω, Iω, phase; energy, λ0=NaN, mode=:lowest, polarisation=:linear, propagator=nothing)')
                Main.duv = Main.eval('duv = prop_capillary(radius, flength, gas, pressure; λ0, pulses=pulse, trange=400e-15, λlims=(150e-9, 4e-6))')
                
        else:
            Main.duv = Main.eval('duv = prop_capillary(radius, flength, gas, pressure; λ0, τfwhm, energy, trange=400e-15, λlims=(150e-9, 4e-6))')

        Main.eval('t, Et = Processing.getEt(duv)')
        Main.eval("λ, Iλ = Processing.getIω(duv, :λ, flength)")
        λ = Main.λ
        Iλ = Main.Iλ
        t = Main.t
        Et_allz=Main.Et #array of Et at all z 
        Et=Et_allz[:,-1] #last item in each element is pulse shape at the end
        Et0=Et_allz[:,0]

        plt.figure()
        plt.plot(λ*10**9,Iλ)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Spectral Energy Density (J/m)")
        plt.title('Output Spectrum from the Fibre produced by Optimal Parameters')

        plt.figure()
        plt.plot(t,Et,label="Output Pulse")
        plt.plot(t,Et0,label="Input Pulse")
        plt.xlabel("Time (s)")
        plt.ylabel("Electric Field (a.u.)")
        plt.title('Input and Output Pulses in the Time Domain')
        plt.legend()
        plt.show()

    return optimizer.max, optimizer.res
