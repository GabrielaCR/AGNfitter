"""

%%%%%%%%%%%%%%%%%%

MCMC_AGNfitter.py

%%%%%%%%%%%%%%%%%%

This script contains 

"""

import emcee #Author: Dan Foreman-Mackey (danfm@nyu.edu)
import sys,os
import time
import pickle
from . import PARAMETERSPACE_AGNfitter as parspace
import ultranest
from ultranest import ReactiveNestedSampler, stepsampler, dychmc, popstepsampler
import numpy as np
from ultranest.plot import cornerplot, traceplot


if __name__ == 'main':
    main(sys.argv[1:])



"""==================================================
 SAMPLING FUNCTIONS
=================================================="""

def main(data, models, P, mc):

    """
    Main function for the MCMC sampling.
    Integrates Emcee with parameter settings and mcmc settings.
    Saves chains into files.

    ##input:
    - object data of class DATA (DATA_AGNfitter.py)
    - dictionary P, of parameter settings (PARAMETERSPACE_AGNfitter.py)
    - dictionary mc, of mcmc settings (RUN_AGNfitter_multi.py)
    """

    print( '......................................................')
    print( 'model parameters', P['names'])
    print( 'minimum values', ['{0:.2f}'.format(k) for k in P['min']])
    print( 'maximum values', ['{0:.2f}'.format(k) for k in P['max']])
    print( '......................................................')
    print ( mc['Nwalkers'], 'walkers')

    Npar = len(P['names'])
    
    # Needed for ultranest
    #print(parspace.ymodel(data.nus, data.z, data.dlum, models, P, *pars))


    if not os.path.lexists(data.output_folder+str(data.name)):
        os.mkdir(data.output_folder+str(data.name))

    def my_posterior(params):
        posterior = []

        for i in range(params.shape[0]):  #For each test sample
            posterior.append(parspace.ln_probab(tuple(params[i]), data, models, P))

        return np.array(posterior)


    def my_params_transform(cube):
        params = cube.copy()

        for i in range(params.shape[0]):  #For each test sample
            normalization = []
            for k in range(params.shape[1]):
                normalization.append( (P['max'][k]-P['min'][k])*cube[i][k] + P['min'][k])
            params[i] = normalization
        return params

    # Ultranest
    sampler = ultranest.ReactiveNestedSampler(P['names'], my_posterior, my_params_transform, resume=True, log_dir=  data.output_folder+str(data.name)+'/ultranest', vectorized=True)

    sampler.stepsampler = ultranest.stepsampler.SliceSampler( nsteps=10,
         generate_direction=ultranest.stepsampler.generate_mixture_random_direction)  # step sampling technique for high dimensional spaces
    #without this, my code can take more than 2 hours for 3 sources

    #sampler.stepsampler = ultranest.dychmc.DynamicCHMCSampler(scale = 1, nsteps = 20, adaptive_nsteps=False, delta=0.9, nudge=1.04)
    #sampler.stepsampler =  ultranest.popstepsampler.PopulationSliceSampler(popsize = 10, nsteps = 20, generate_direction = ultranest.stepsampler.generate_mixture_random_direction, scale=1.0, scale_adapt_factor=0.9, log=False, logfile=None)

    sampler.run(
     min_num_live_points=400,
     min_ess=400, # number of effective samples
     max_num_improvement_loops=3, # how many times to go back and improve
 ) 

    sampler.plot_run()
    sampler.plot_trace()
    sampler.plot_corner()




