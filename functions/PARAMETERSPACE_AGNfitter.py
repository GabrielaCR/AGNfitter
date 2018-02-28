
"""%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAMETERSPACE_AGNfitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains functions that rule the exploration of the parameter space.

It contains:



"""
from __future__ import division
import pylab as pl
import numpy as np
from math import pi
import time
from collections import Iterable
import itertools
import pickle
import MODEL_AGNfitter as model

def flatten(lis):
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, basestring):
             for x in flatten(item):
                 yield x
         else:        
             yield item

def Pdict (data):

    """Constructs a dictionary P with the name of all model parameters as keys. 
    The values are tuples with the same length, with the parameter limits. 

    ## inputs:
    - data object
    """
    P = dict()
    ga,sb,to,bb= data.dictkey_arrays
    par_mins = list(flatten([min(i.astype(float)) if i.ndim==1 else [min(i[j].astype(float)) for j in range(len(i))] \
               for i in [np.array(ga.pars_modelkeys), np.array(sb.pars_modelkeys), np.array(to.pars_modelkeys), np.array(bb.pars_modelkeys)]]))
    par_maxs = list(flatten([max(i.astype(float)) if i.ndim==1 else [max(i[j].astype(float)) for j in range(len(i))] \
               for i in [np.array(ga.pars_modelkeys), np.array(sb.pars_modelkeys), np.array(to.pars_modelkeys), np.array(bb.pars_modelkeys)]]))
    
    ## Add normalization parameters:
    if len(bb.par_names)==1:
        [par_mins.append(i) for i in [0.,0.,0.,0.]]
        [par_maxs.append(i)for i in [15.,15.,15.,15.]]
        normpars=['GA','SB','TO','BB'] 

    else:
        [par_mins.append(i) for i in [0.,0.,0.]]
        [par_maxs.append(i)for i in [15.,15.,15.]]           
        normpars=['GA','SB','TO'] 

    all_pars = list(itertools.chain.from_iterable([ ga.par_names, sb.par_names,to.par_names, bb.par_names ,normpars]))  
    npc= [len(ga.par_names),len(sb.par_names),len(to.par_names), len(bb.par_names), len(normpars)]  

    P['names'] = all_pars
    Npar = len(P['names'])
    P['priortype'] = np.array(['non-info']*Npar)
    P['min'] = par_mins
    P['max'] = par_maxs    
    P['idxs'] = [0, sum(npc[0:1]),sum(npc[0:2]),sum(npc[0:3]),sum(npc[0:4])]

    return P  

"""-------------------------
PRIOR, LIKELIHOOD, POSTERIOR 
---------------------------"""


def ln_prior(z, dlum,bands, gal_Fnu, P, pars):

    """Calculates the prior probability on the parameters.

    (1) Flat prior on all parameters considerinf limits described in P.
    (2) Flat prior on the galaxy, using the B-band magnitude expected 
         from the galaxy luminosity function as a maximum of the prior.

    ## inputs: z, dlum,bands, gal_Fnu, P, pars
    ## output: -inf or 0.
    """

    #1. Flat priors
    for i,p in enumerate(pars):
        if P['priortype'][i]=='non-info' and not (P['min'][i] < p < P['max'][i]):
            return -np.inf

    #2. Prior on the luminosity
    B_band_expected, B_band_thispoint = galaxy_Lumfct_prior( z, dlum, bands, gal_Fnu)    # Bband expectations
    #if Bband magnitude in this trial is brighter than expected by the luminosity function, dont accept this one
    if B_band_thispoint < (B_band_expected - 5):#2.5):
        return -np.inf

    return 0.



def ln_likelihood(x, y, ysigma, z, ymodel):

    """Calculates the likelihood function.

    It includes the restriction of taking into account only 
    frequencies lower than the Ly-alpha line, to ensure being free
    IGM absorption. (specially relevant for z>3.)

    ## inputs:
    - x, y, ysigma, z
    - y_model calculate with fct ymodel()

    ## output:
    - (-1 * ln(likelihood))"""
    #x_valid:
    #only frequencies with existing data (no detections nor limits F = -99)        
    #Consider only data free of IGM absorption. Lyz = 15.38 restframe        
    x_valid = np.arange(len(x))[(x< np.log10(10**(15.38)/(1+z))) & (y>-99.)]

    resid = [(y[i] - ymodel[i])/ysigma[i] for i in x_valid]
    return -0.5 * np.dot(resid, resid)




def ln_probab(pars, data, P):

    """Calculates the posterior probability as Ppos= Pprior + Pdata
    ## inputs:
    - pars
    - object data
    - dictionary P

    ## output:
    - POSTERIOR probabiliy

    ## dependencies:
    - MCMC_AGNfitter.py"""

    y_model, bands, gal_Fnu = ymodel(data.nus,data.z, data.dlum, data.dictkey_arrays,data.dict_modelfluxes, P, *pars)

    lnp = ln_prior(data.z, data.dlum, bands,gal_Fnu, P, pars)

    if np.isfinite(lnp):    
        posterior = lnp + ln_likelihood(data.nus,data.fluxes,data.fluxerrs, data.z, y_model)     
        return posterior
    return -np.inf



"""------------------------------------
CONSTRUCT TOTAL MODEL 
------------------------------------"""

def ymodel(data_nus, z, dlum, dictkey_arrays, dict_modelfluxes, P, *par):

    """Constructs the total model from parameter values.

    ## inputs: data_nus, z, dictkey_arrays, dict_modelfluxes, *par

    ## output:
    - total model
    - bands
    - galaxy_flux (to be used by ln_prior)

    ## dependencies:
    This fct is used 
    (1)in ln_likelihood, this same script.
    (2)in scriptPLOTandWRITE.

    """
    t0 = time.time()

    STARBURSTFdict , BBBFdict, GALAXYFdict, TORUSFdict,_,_,_,_,_,_= dict_modelfluxes

    gal_obj,sb_obj,tor_obj, bbb_obj = dictkey_arrays
 
    par =par[0:len(par)]

    ## Use pick_nD if model has more than one parameter,
    ## and pick_1D if it has only one.
    gal_obj.pick_nD(par[P['idxs'][0]:P['idxs'][1]])  
    sb_obj.pick_nD(par[P['idxs'][1]:P['idxs'][2]]) 
    tor_obj.pick_1D(par[P['idxs'][2]:P['idxs'][3]])            

    ### Use when BBB model has only 1 par (Richards 06)
    # GA, SB, TO, BB = par[-4:]
    # bbb_obj.pick_1D(par[P['idxs'][3]:P['idxs'][4]])

    ### Use when the BBB model has >1 parameters
    GA, SB, TO = par[-3:]
    BB = 0.
    bbb_obj.pick_nD(par[P['idxs'][3]:P['idxs'][4]])

    t2 = time.time()

    try: 
        bands, gal_Fnu = GALAXYFdict[tuple(gal_obj.matched_parkeys)]   
        _, sb_Fnu= STARBURSTFdict[tuple(sb_obj.matched_parkeys)] 
        _, bbb_Fnu = BBBFdict[tuple(bbb_obj.matched_parkeys)]   ### use when the BBB model has >1 parameters
        #_, bbb_Fnu = BBBFdict[bbb_obj.matched_parkeys]   ### use when BBB model has only 1 par (Richards 06)
        _, tor_Fnu= TORUSFdict[tor_obj.matched_parkeys] 
    except ValueError:
        print 'Error: Dictionary does not contain some values'
    t3 = time.time()

    # Renormalize to have similar amplitudes. Keep these fixed!
    sb_Fnu_norm =sb_Fnu/1e20    
    bbb_Fnu_norm = bbb_Fnu/ (dlum)**2
    gal_Fnu_norm = gal_Fnu/1e18
    tor_Fnu_norm = tor_Fnu/ 1e-40

    # Total SED sum
    #--------------------------------------------------------------------

    lum = 10**(SB)* sb_Fnu_norm  + 10**(BB)*bbb_Fnu_norm    \
          + 10**(GA)*gal_Fnu_norm  +(10**TO) *tor_Fnu_norm

    #--------------------------------------------------------------------

    lum = lum.reshape((np.size(lum),))    
    return lum, bands, 10**(GA)*gal_Fnu_norm





def galaxy_Lumfct_prior( z, dlum, bands, gal_flux):

    """This function calculates 
    (1)the Bband magnitude of the galaxy template
    (2)the Bband magnitude expected from the galaxy luminosity function
         given in Iovino et al. (2010)
    ## inputs:
    -    z(float), dlum(float), bands(array), gal_flux (array)"""

    # Calculated B-band at this parameter space point
    h_70 = 1.
    lumfactor = (4. * pi * dlum**2.)

    flux_B = gal_flux[(14.790 < bands)&(bands < 14.870)]
    if len(flux_B)>1:
        flux_B = flux_B[0]
    mag1= -2.5 * np.log10(flux_B) - 48.6
    distmod = -5.0 * np.log10((dlum/3.08567758e24 *1e6)/10) 
    abs_mag1 = mag1 + distmod
    thispoint1 = abs_mag1

    lum_B = lumfactor * flux_B
    abs_mag = 51.6 - 2.5 *np.log10(lum_B)
    thispoint = abs_mag

    # Expected B-band calculation (Iovino et al. (2010))
    expected = -20.3 - (5 * np.log10(h_70) )- (1.1 * z) 


    return expected,thispoint


"""--------------------------------------
Obtain initial positions
--------------------------------------"""


def get_initial_positions(nwalkers, P):

    """Returns the initial positions.
    ## inputs:
    - number of walkers (int)
    - dictionary P
    ## output:
    - list of positions of length len(P.keys)"""
    Npar = len(P['names']) 
    p0 = np.random.uniform(size=(nwalkers, Npar))

    for i in range(Npar):
        p0[:, i] = 0.5*(P['min'][i] + P['max'][i]) + (2* p0[:, i] - 1) * (1)
    p0[:,8]= 0.1 + (2* p0[:, 8] - 1) * (0.001)
    
    return p0


def get_best_position(filename, nwalkers, P):

    """Returns the best positions after burn-in phases.
    ## inputs:
    - filename (str), nwalkers(int), P (dict)
    ## output:
    - list of positions of length (P.keys)"""

    Npar = len(P['names']) 
    #all saved vectors    
    f = open(filename, 'rb')
    samples = pickle.load(f)
    f.close()

    #index for the largest likelihood     
    i = samples['lnprob'].ravel().argmax()
    #the values for the parameters at this index
    P['ml'] = samples['chain'].reshape(-1, Npar)[i]

    p1 = np.random.normal(size=(nwalkers, Npar))

    for i in range(Npar):
        p =P['ml'][i]
        
        p1[:, i] = p + 0.00001 * p1[:, i]

    return p1




