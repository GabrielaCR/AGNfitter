"""

%%%%%%%%%%%%%%%%%%

MCMC_AGNfitter.py

%%%%%%%%%%%%%%%%%%

This script contains 

"""

import emcee #Author: Dan Foreman-Mackey (danfm@nyu.edu)
import sys,os
import time
import numpy as np
import cPickle
import PARAMETERSPACE_AGNfitter as parspace
from DATA_AGNfitter import DATA



def main(data, P, mc):

    """
    Main function for the MCMC sampling.
    Integrates Emcee with parameter settings and mcmc settings.
    Saves chains into files.

    ##input:
    - object data of class DATA (DATA_AGNfitter.py)
    - dictionary P, of parameter settings (PARAMETERSPACE_AGNfitter.py)
    - dictionary mc, of mcmc settings (RUN_AGNfitter_multi.py)
    """

    path = os.path.abspath(__file__).rsplit('/', 1)[0]

    print '......................................................'
    print 'model parameters', P['names']
    print 'minimum values', ['{0:.2f}'.format(k) for k in P['min']]
    print 'maximum values', ['{0:.2f}'.format(k) for k in P['max']]
    print '......................................................'
    print mc['Nwalkers'], 'walkers'

    Npar = len(P['names'])


    sampler = emcee.EnsembleSampler(
            mc['Nwalkers'], Npar, parspace.ln_probab,
            args=[data, P],  daemon= True)


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
        print '%.2g min elapsed' % ((time.time() - t1)/60.)


    ## MCMC SAMPLING ##
    if mc['Nmcmc'] > 0:

        t2 = time.time()
        run_mcmc(sampler, p_maxlike, data.name,data.output_folder, mc)
        print '%.2g min elapsed' % ((time.time() - t2)/60.)
    del sampler.pool    


if __name__ == 'main':
    main(sys.argv[1:])



"""==================================================
 SAMPLING FUNCTIONS
=================================================="""


def run_burn_in(sampler, mc, p0, sourcename, folder, setnr):
    """ Run and save a set of burn-in iterations."""

    print 'Running burn-in nr. '+ str(setnr)+' with %i steps' % mc['Nburn']
    
    iprint = mc['iprint']

    # note the results are saved in the sampler object.
    for i,(pos, lnprob, state) in enumerate(sampler.sample(p0, iterations=mc['Nburn'])):
        i += 1
        if not i % iprint:
            print i
        
    save_chains(folder+str(sourcename)+'/samples_burn1-2-3.sav', sampler, pos, state)

    return pos, state   


def run_mcmc(sampler, pburn, sourcename, folder, mc):
    """
    Run MCMC sampling and save.    
    """

    sampler.reset()

    iprint = mc['iprint']
    print "Running MCMC with %i steps" % mc['Nmcmc']

    for i,(pos, lnprob, state) in enumerate(sampler.sample(pburn, iterations=mc['Nmcmc'])): 
        i += 1
        if not i % iprint:
            print i
            
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
    cPickle.dump(dict(
        chain=sampler.chain, accept=sampler.acceptance_fraction,
        lnprob=sampler.lnprobability, final_pos=pos, state=state, acor=sampler.acor), f, protocol=2)
    f.close()





