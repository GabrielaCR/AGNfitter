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
from ultranest.stepsampler import OrthogonalDirectionGenerator
import importlib


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

    # Functions for ultranest
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

    if mc['sampling_algorithm'] == 'ultranest':

        if not os.path.lexists(data.output_folder+str(data.name)):
            os.mkdir(data.output_folder+str(data.name))

        sampler = ultranest.ReactiveNestedSampler(P['names'], my_posterior, my_params_transform, resume=True, log_dir=  data.output_folder+str(data.name)+'/ultranest', vectorized=True)

        if mc['direction_generation'] == 'de-mix':
            direction_stepsampler = ultranest.stepsampler.generate_mixture_random_direction
        elif mc['direction_generation'] == 'region-slice':
            direction_stepsampler = ultranest.stepsampler.generate_region_oriented_direction
        elif mc['direction_generation'] == 'cube-ortho-harm':
            direction_stepsampler = OrthogonalDirectionGenerator(ultranest.stepsampler.generate_random_direction)


        sampler.stepsampler = ultranest.stepsampler.SliceSampler( nsteps=20,
         generate_direction = direction_stepsampler)  # step sampling technique for high dimensional spaces
                                                      #without this, the code can take more than 2h for 3 sources

        sampler.run( min_num_live_points= mc['live_points'], min_ess= mc['min_ess'], max_num_improvement_loops = mc['num_loops'] , ) 
        sampler.plot_run()
        sampler.plot_trace()
        sampler.plot_corner()

    elif mc['sampling_algorithm'] == 'emcee':
        #Change default value of quiet = False in emcee auto correlation time function
        file_autocorr = os.path.dirname(emcee.__file__) + '/autocorr.py'

        acor_2r = open(file_autocorr, 'r')
        Lines = acor_2r.readlines()
        acor_2r.close()
        with open(file_autocorr, 'w') as acor_py:
            for line in Lines:
                if 'integrated_time(x, c=5, tol=50, quiet=False)' in line:
                    line = line.replace('quiet=False', 'quiet=True')
                    print('Autocorr python file has been changed. Quiet input of integrated_time became True')
                acor_py.write(line)
        acor_py.close()
        importlib.reload(emcee.autocorr)

        sampler = emcee.EnsembleSampler( mc['Nwalkers'], Npar, parspace.ln_probab, args=[data, models, P])

        ## BURN-IN SETS ##
        if mc['Nburn'] > 0:
            t1 = time.time()
            if not os.path.lexists(data.output_folder+str(data.name)):
                os.mkdir(data.output_folder+str(data.name))

            p_maxlike = parspace.get_initial_positions(mc['Nwalkers'], P)
            Nr_BurnIns = mc['Nburnsets']  

            for i in range(Nr_BurnIns):
                p_maxlike, state = run_burn_in(sampler, mc, p_maxlike, data.name, data.output_folder, i)
                savedfile = data.output_folder+str(data.name)+'/samples_burn1-2-3.sav'
                p_maxlike = parspace.get_best_position(savedfile, mc['Nwalkers'], P)
            print( '%.2g min elapsed' % ((time.time() - t1)/60.))

        ## MCMC SAMPLING ##
        if mc['Nmcmc'] > 0:
            t2 = time.time()
            run_mcmc(sampler, p_maxlike, data.name,data.output_folder, mc)
            print( '%.2g min elapsed' % ((time.time() - t2)/60.))
        del sampler.pool  

        #Change to default value of quiet = False in emcee auto correlation time function

        acor_2r = open(file_autocorr, 'r')
        Lines = acor_2r.readlines()
        acor_2r.close()
        with open(file_autocorr, 'w') as acor_py:
            for line in Lines:
                if 'integrated_time(x, c=5, tol=50, quiet=True)' in line:
                    line = line.replace('quiet=True', 'quiet=False')
                    print('Autocorr python file has been changed. Quiet input of integrated_time became False again')
                acor_py.write(line)
        acor_py.close()

"""==================================================
 SAMPLING FUNCTIONS
=================================================="""


def run_burn_in(sampler, mc, p0, sourcename, folder, setnr):
    """ Run and save a set of burn-in iterations."""

    print( 'Running burn-in nr. '+ str(setnr)+' with %i steps' % mc['Nburn'])
    
    iprint = mc['iprint']

    # note the results are saved in the sampler object.
    for i,(pos, lnprob, state) in enumerate(sampler.sample(p0, iterations=mc['Nburn'])):
        i += 1
        if not i % iprint:
            print( i )
        
    save_chains(folder+str(sourcename)+'/samples_burn1-2-3.sav', sampler, pos, state)

    return pos, state   


def run_mcmc(sampler, pburn, sourcename, folder, mc):
    """
    Run MCMC sampling and save.    
    """

    sampler.reset()

    iprint = mc['iprint']
    print( "Running MCMC with %i steps" % mc['Nmcmc'])

    for i,(pos, lnprob, state) in enumerate(sampler.sample(pburn, iterations=mc['Nmcmc'])): 
        i += 1
        if not i % iprint:
            print( i)
            
    save_chains(folder+str(sourcename)+'/samples_mcmc.sav', sampler, pos, state)   


def save_chains(filename, sampler, pos, state):
    """
    Save dictionary which contains:
    -chains
    -acceptance_fraction
    -lnprob
    -last positions
    -autocorrelation time 
    into .sav files, using cPickle.
    """
    f = open(filename, 'wb')
    pickle.dump(dict(
        chain=sampler.chain, accept=sampler.acceptance_fraction,
        lnprob=sampler.lnprobability, final_pos=pos, state=state, acor=sampler.acor), f, protocol=2)
    f.close()

