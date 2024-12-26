# General imports
import os
import numpy as np
from astropy.cosmology import FlatLambdaCDM

# Setting cosmology
cosmo = FlatLambdaCDM(H0=73, Om0=0.3, Tcmb0=2.725)

# Xspec, BXA, and UltraNest imports and config
from xspec import *
Xset.chatter = 1  # Controls pyXPSEC output
Xset.cosmo = "73"  # Cosmology for pyXPSEC
import bxa.xspec as bxa
bxa.BXASolver.allowed_stats.append("chi")  # Enable Gaussian statistics
from bxa.xspec.solver import set_parameters  # For later use
from ultranest.plot import PredictionBand
Fit.statMethod = "chi"  # Set pyXPSEC statistics to Gaussian
Fit.query = "yes"
Plot.device = "/xs"  # Plotting console for pyXPSEC
Plot.yLog = True
Plot.background = True

# Load diskSED model
from diskSED import diskSED_info, diskSED
AllModels.addPyMod(diskSED, diskSED_info(), 'add')  # Additive model

# Load reddenSF model (Calzetti et al. attenuation law)
from diskSED import reddenSF_info, reddenSF
AllModels.addPyMod(reddenSF, reddenSF_info(), 'mul')  # Multiplicative model

# Load reddenCCM model (CCM extinction law) This is the same as the XSPEC native 'redden' but extended to the Lyman limit.
from diskSED import reddenCCM_info, reddenCCM
AllModels.addPyMod(reddenCCM, reddenCCM_info(), 'mul')  # Multiplicative model

# Current working directory
cwd = os.getcwd()

# Source redshift
z = 0.0206

# Galactic column density at 14li's position (units of 10^22)
N_H_G = 1.95E+20 / 1e22

# Quick cleaning (in case something has been loaded already)
AllData.clear()
AllModels.clear()
AllData.clear()

Plot.xAxis = "keV"
# Set energy bins for the model (working from 5e-4 keV to 2 keV with 400 log bins)
AllModels.setEnergies('0.0005 2.0 400 log')

# Load X-ray and UV/optical data
E_lim = [0.25, 0.9]  # Energy limits for X-ray spectrum, use wathever criteria you think is resonable
os.chdir(os.path.join(cwd, 'data', 'SED'))
x_ray_file = '14li_E1_grp.fits'
uvopt_file = 'uvopt.pha'
AllData('1:1 ' + str(x_ray_file) + ' 1:2 ' + str(uvopt_file))

# Ignore bad channels
Xray = AllData(1)
Xray.response.arf = '14li_E1_arf.fits' # Sometimes we need to add the arf manually, it fails to find automatically.
Xray.ignore('0.0-' + str(E_lim[0]) + ' ' + str(E_lim[1]) + '-**')
uvopt = AllData(2)
uvopt.response = 'uvopt.rmf'
AllData.ignore("bad")

#If you had another epoch you could use this to load it:
'''
E_lim  = [0.25, .8]
os.chdir(os.path.join(cwd, 'data', 'SED2'))
x_ray_file2 = '14li_E2_grp.fits'
uvopt_file2 = 'uvopt2.pha'
AllData('2:3 '+str(x_ray_file2) + ' 2:4 ' +str(uvopt_file2))
Xray2 =  AllData(3)
Xray2.ignore('0.0-'+ str(E_lim[0])+ ' '+ str(E_lim[1])+ '-**')
uvopt2 = AllData(4)
uvopt2.response = 'uvopt2.rmf'
'''

# Let's Plot the data
Plot("data") # Yes, the units look weird this way. But is just to make sure, it is all loaded. You see an X-ray spectrum (with background) and an opt/UV SED.

os.chdir(cwd)

# Set the model

m = Model('phabs*redden*zashift*phabs*reddenSF*(diskSED)') 
# Please read arXiv:2408.17296 for a detailed description:
# phabs*reddenCCM accounts for Galactic X-ray neutral gas absorption and Galactic dust extinction, respectively.
# zashift redshifts the model.
# phabs*reddenSF accounts for intrinsic X-ray neutral gas absorption and intrinsic dust attenuation, respectively.
# Note: We use phabs instead of, e.g., tbabs, because tbabs corrects down to UV wavelengths (not only X-ray) and using it would result in double correction in the UV.

AllModels.show() # The model components are shown below.

# Let's edit these parameters, fix what needs to be fixed, and add reasonable ranges to those that will be fitted.

# phabs Nh
AllModels(1)(1).values = (N_H_G, -1)  # Fixing

# reddenCCM
AllModels(1)(2).link =  '1 * 1e22 / 6.88e21' # This assumes a Galactic Gas-to-dust Ratio with R_v  = 3.1. See comments below.

# zashift Redshift
AllModels(1)(3).values = (z, -1)  # Fixing

# phabs Nh
AllModels(1)(4).values = (0.04, 0.001, 0.02, 0.02, 0.1, 0.1)  
# Check the pyxspec manual for the meaning of each of the six entries. The first entry is the initial guess, but that is not used in a Bayesian fitting.

# redden
AllModels(1)(5).link = '4 * 1e22 / 8.9505e21'  
# Linking intrinsic E(B-V) (parameter 5) to intrinsic N_H (parameter 4). 
# Assumes a Galactic Gas-to-Dust Ratio (https://arxiv.org/abs/0903.2057) with R_v = 4.05 (consistent with Calzetti's law).
# You could also make this free or use, e.g., a Gaussian prior centered around an independent estimate 
# for the intrinsic E(B-V), such as from the Balmer decrement of the nuclear optical spectrum of the source.

# Parameters for diskSED. Feel free to modify give the expectations from your source/observations; note that larger ranges will increase convergence time.

# diskSED R_in*: Expected range for the Rin* parameter in km.
AllModels(1)(6).values = (3e7, 10, 3e7, 3e7, 1e8, 1e8)

# diskSED T_p: Expected range for inner disk temperature in Kelvin.
AllModels(1)(7).values = (2e5, 1e-2, 1e5, 1e5, 1e6, 1e6)

# diskSED R_ratio: Expected outer-to-inner radius ratio range.
AllModels(1)(8).values = (10, 1e-1, 5, 5, 30, 30)

# diskSED Distance: Distance to the source in Mpc. Needs to be fixed.
AllModels(1)(9).values = (cosmo.luminosity_distance(z).to('Mpc').value, -1)

# diskSED norm: This must always be set to 1.
AllModels(1)(10).values = (1, -1)  
# Rin* is already the normalization parameter of the model.

# Check the model below; there should be exactly 4 free parameters, and 2 Tied.
AllModels.show()

# Using BXA:

# These creates the priors
prior_NH = bxa.create_uniform_prior_for(AllModels(1), AllModels(1)(4)) # The number 1 is the data group (here we only have 1), while the other number is the paramter order, e.g. 4, see above. See BXA tutorial for details.
prior_R = bxa.create_loguniform_prior_for(AllModels(1), AllModels(1)(6))
prior_T = bxa.create_uniform_prior_for(AllModels(1), AllModels(1)(7))
prior_R_ratio = bxa.create_uniform_prior_for(AllModels(1), AllModels(1)(8))

# This creates all the rest
solver = bxa.BXASolver(transformations=[prior_NH, prior_R, prior_T, prior_R_ratio])

# Define/load here
paramnames_bxa = solver.paramnames
prior_bxa  = solver.prior_function
loglike_bxa  = solver.log_likelihood
prior_transform_bxa = solver.transformations

print('------------------------------------------')
print(paramnames_bxa) # This is just a list with name of paramters

# OK, now we can run some Baeysian inference. I will use UltraNest.
import ultranest
from ultranest import ReactiveNestedSampler
import pandas as pd
# We wanna make pyXSPEC quiet now. So just Ultranest can display progress.
Xset.chatter = 1

os.chdir(os.path.join(cwd, 'data')) # make sure we are where we want.
outputfiles_basename = 'diskSED_run' # Name of the folder we gonna save the sampling results

#picking BXA funcitons
sampler = ReactiveNestedSampler(paramnames_bxa, loglike_bxa, prior_bxa, log_dir=outputfiles_basename, resume='overwrite')

#picking our own funcitons
#sampler = ReactiveNestedSampler(paramnames, loglike, prior, log_dir=outputfiles_basename, resume='overwrite')

result = sampler.run(frac_remain=1e-2, min_num_live_points=200,min_ess=400) # See UltraNest for details on what the varaibles mean. These are standard values. 
# Some intuition: 
# decreasing frac_remain (e.g. to 5e-2, or 1e-1) wil make the sampling finish earlier, running longer does not mean necessarily a better fit.
# decreasing min_num_live_points decreases the 'resolution' of the sampling, but makes the computations faster.
sampler.print_results()


posterior_df = pd.DataFrame(data=result["samples"], columns=solver.paramnames)
posterior_df.to_csv(os.path.join(cwd, 'data', "%s/posterior.csv" %(outputfiles_basename)), index=False) # Just saving the posterior into a .csv file, cause it is easier.

# This will run and save the results inside a 'diskSED_run' folder.
# It will take time (10-20 min). If you wanna speed up things, e.g. by running in multiple cores. There will be a 'fitting.py' script, which we can parallelize with MPI (https://www.open-mpi.org).
# In the terminal, inside the examples folder, just do 'mpiexec -N X python3 fitting.py' where X is the number of cores you wanna run the sampling.