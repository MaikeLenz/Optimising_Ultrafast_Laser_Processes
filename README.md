# Optimising_Ultrafast_Laser_Processes
Final Code, MSci Project, Molly Smith and Maike Lenz

This code can be used to optimise hollow-core fibre experiments to produce ultrashort laser pulses. 

It relies on the installation of two packages:
1) https://github.com/LupoLab/Luna
The Luna package simulates the physics that occurs within a hollow-core fibre.
2) https://github.com/fmfn/BayesianOptimization
The Bayesian optimisation module is used to optimise the experiment.
For details on the installation of these packages and our notes on the installation process please see the getting started section below and the respective github repositries.

## Getting Started
### To install Luna:
1) Download Julia https://julialang.org/downloads/
2) In cmd, run:
	* julia
	* ]
	* add https://github.com//CoolProp/CoolProp.jl 
	* add https://github.com/LupoLab/Luna
	* Ctrl c
	* using Luna
	* output =prop_capillary(125e-6, 3, :He, 1; λ0=800e-9, energy=120e-6, τfwhm=10e-15, λlims=(150e-9, 4e-6), trange=1e-12)
	* output["Eω"]

	This should precompile Luna and run the simulation through julia - for more info https://github.com/chrisbrahms/Luna.jl. Precompiling Luna may take a while.
Luna documentation: http://lupo-lab.com/Luna.jl/dev/model/model.html
3) I had some issues precompiling Luna due to the julia module PyCall not being setup right. Run in cmd:
	* julia
	* ENV["PYTHON"]="Python_Path_Name"
	* ]
	* build PyCall
	Then try adding Luna again
4) In python: pip install julia
5) In python, run:
	* import julia
	* julia.install()
	Need to be careful that the python version used is the same one that has been set in julia through ENV["PYTHON"]
6) Run script Sim_python.py through cmd to run julia simulation through python

### To install Bayesian Optimisation:
1) Clone https://github.com/fmfn/BayesianOptimization 
2) Run the file setup.py install
3) pip install bayesian-optimization
4) In python, run:
	* from bayes_opt import BayesianOptimization

If this doesn't run, use sys to add the path to the Bayesian Optimization module to your python code. Then check that the module bayesian-optimization is installed in the version of python that you are using.
New issue: BayesianOptimization currently only works with scipy=1.7 or earlier.
