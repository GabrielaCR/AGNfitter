
"""%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAMETERSPACE_AGNfitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains functions that rule the exploration of the parameter space.

It contains:



"""
from __future__ import division
import numpy as np
import time
from collections import Iterable
import itertools
import pickle
from . import PRIORS_AGNfitter as priors

def flatten(lis):
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, str):
             for x in flatten(item):
                 yield x
         else:        
             yield item

def Pdict (data, models):

    """Constructs a dictionary P with the name of all model parameters as keys. 
    The values are tuples with the same length, with the parameter limits. 

    ## inputs:
    - data object
    """
    P = dict()
    ga,sb,to,bb, agnrad= models.dictkey_arrays
    if models.settings['RADIO'] == True and (agnrad.pars_modelkeys != ['-99.9']).all() : #If there is a radio model with fitting parameters
        par_mins = list(flatten([min(i.astype(float)) if i.ndim==1 else [min(i[j].astype(float)) for j in range(len(i))] \
               for i in [np.array(ga.pars_modelkeys), np.array(sb.pars_modelkeys), np.array(to.pars_modelkeys), np.array(bb.pars_modelkeys), np.array(agnrad.pars_modelkeys)]]))
        par_maxs = list(flatten([max(i.astype(float)) if i.ndim==1 else [max(i[j].astype(float)) for j in range(len(i))] \
               for i in [np.array(ga.pars_modelkeys), np.array(sb.pars_modelkeys), np.array(to.pars_modelkeys), np.array(bb.pars_modelkeys), np.array(agnrad.pars_modelkeys)]]))

    else:            #If there isn't radio data available or if the model have fix parameters, there aren't included in the exploration of parameters space
        par_mins = list(flatten([min(i.astype(float)) if i.ndim==1 else [min(i[j].astype(float)) for j in range(len(i))] \
               for i in [np.array(ga.pars_modelkeys), np.array(sb.pars_modelkeys), np.array(to.pars_modelkeys), np.array(bb.pars_modelkeys)]]))
        par_maxs = list(flatten([max(i.astype(float)) if i.ndim==1 else [max(i[j].astype(float)) for j in range(len(i))] \
               for i in [np.array(ga.pars_modelkeys), np.array(sb.pars_modelkeys), np.array(to.pars_modelkeys), np.array(bb.pars_modelkeys)]]))

    ## Add normalization parameters:
    [par_mins.append(i) for i in [-10,-10.,-10]]
    [par_maxs.append(i)for i in [ 10, 10,10]]  ##ALTERADO
    normpars=['GA','SB','TO'] 

    if models.settings['BBB']=='R06': #Only the Richards et al. model need normalization parameter
        par_mins.append(-10.)
        # par_maxs.append(10)  
        par_maxs.append(-9)  
        normpars.append('BB') 

    if models.settings['RADIO'] == True:
        par_mins.append(-10.)
        par_maxs.append(10)           
        normpars.append('RAD')

    if models.settings['RADIO'] == True and (agnrad.pars_modelkeys != ['-99.9']).all(): #If there is a radio model with fitting parameters
        all_pars = list(itertools.chain.from_iterable([ ga.par_names, sb.par_names,to.par_names, bb.par_names, agnrad.par_names, normpars]))  
        npc= [len(ga.par_names),len(sb.par_names),len(to.par_names), len(bb.par_names), len(agnrad.par_names), len(normpars)]  
        P['priortype'] = [ga.par_types,sb.par_types, to.par_types, bb.par_types, agnrad.par_types,['free']*len(normpars)] 
        P['idxs'] = [0, sum(npc[0:1]),sum(npc[0:2]),sum(npc[0:3]),sum(npc[0:4]), sum(npc[0:5])]

    else:           #If there isn't radio data available or if the model have fix parameters, there aren't included in the exploration of parameters space
        all_pars = list(itertools.chain.from_iterable([ ga.par_names, sb.par_names,to.par_names, bb.par_names ,normpars]))  
        npc= [len(ga.par_names),len(sb.par_names),len(to.par_names), len(bb.par_names), len(normpars)]  
        P['priortype'] = [ga.par_types,sb.par_types, to.par_types, bb.par_types, ['free']*len(normpars)] 
        P['idxs'] = [0, sum(npc[0:1]),sum(npc[0:2]),sum(npc[0:3]),sum(npc[0:4])]

    P['names'] = all_pars
    P['min'] = par_mins
    P['max'] = par_maxs    
    for i,ps in enumerate(P['names']): ### Maximum age is the age of the Universe
        if P['names'][i]=='age':
            P['max'][i] =np.log10(priors.maximal_age(data.z))
        if P['names'][i]=='tau':
            P['max'][i] =np.log10(priors.maximal_age(data.z))
        if P['names'][i]=='Tdust':
            P['max'][i] = 42.
        if P['names'][i]=='EBVgal':####
           P['min'][i] =0.05 ### can be generally assumed for galaxies with log M*> 9.5

    return P  

"""-------------------------
PRIOR, LIKELIHOOD, POSTERIOR 
---------------------------"""


def ln_prior(data, models, P, *pars):

    """Calculates the prior probability on the parameters.

    (1) Flat prior on all parameters considerinf limits described in P.
    (2) Flat prior on the galaxy, using the B-band magnitude expected 
         from the galaxy luminosity function as a maximum of the prior.

    """
    for i,p in enumerate(pars):
        if not (P['min'][i] < p < P['max'][i]):
            return -np.inf

    prior= priors.PRIORS(data, models, P, *pars)
    return prior


def ln_likelihood(x, y, ysigma, z, ymodel, models):

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
    if models.settings['XRAYS'] == 'Prior':                                  #Ignore X-rays data in the likelihood (already taken into account in prior)      
        x_valid = np.arange(len(x))[(x< np.log10(10**(15.38)/(1+z))) & (y>-99.e-23)]
    else:                                                                    #Ignore only UV data because of IGM absorption  
        x_valid = np.arange(len(x))[(x< np.log10(10**(15.38)/(1+z))) | (x > np.log10(10**(16.685)/(1+z))) & (y>-99.e-23)]
    resid = [(y[i] - ymodel[i])/ysigma[i] for i in x_valid]
    return -0.5 * np.dot(resid, resid)


def ln_probab(pars, data, models, P):

    """Calculates the posterior probability as Ppos= Pprior + Pdata
    ## inputs:
    - pars
    - object data
    - dictionary P

    ## output:
    - POSTERIOR probabiliy

    ## dependencies:
    - MCMC_AGNfitter.py"""

    y_model, bands  = ymodel(data.nus, data.z, data.dlum, models, P, *pars)
    lnp = ln_prior(data, models, P, *pars)

    if np.isfinite(lnp): 
        posterior = lnp + ln_likelihood(data.nus, data.fluxes, data.fluxerrs, data.z, y_model, models) 
        return posterior
    return -np.inf


"""------------------------------------
CONSTRUCT TOTAL MODEL 
------------------------------------"""

def ymodel(data_nus, z, dlum, models, P, *par):

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

    par = par[0:len(par)]
    gal_obj,sb_obj,tor_obj, bbb_obj, agnrad_obj = models.dictkey_arrays

    if models.settings['RADIO'] == True:
        if models.settings['BBB']=='R06':
            GA, SB, TO, BB, RAD = par[-5:]
        else:
            GA, SB, TO, RAD = par[-4:]
        if (agnrad_obj.pars_modelkeys != ['-99.9']).all() :             #If there is a radio model with fitting parameters
            agnrad_obj.pick_nD(par[P['idxs'][4]:P['idxs'][5]])
            _, agnrad_Fnu= agnrad_obj.get_fluxes(agnrad_obj.matched_parkeys)
        else:           #If the model have fix parameters, there aren't included in the exploration of parameters space and there is an unique SED template
            all_agnrad_nus, agnrad_Fnu = agnrad_obj.get_fluxes('-99.9')
    else:
        if models.settings['BBB']=='R06':
            GA, SB, TO, BB = par[-4:]
        else:
            GA, SB, TO = par[-3:]

    gal_obj.pick_nD(par[P['idxs'][0]:P['idxs'][1]])  
    sb_obj.pick_nD(par[P['idxs'][1]:P['idxs'][2]]) 
    tor_obj.pick_nD(par[P['idxs'][2]:P['idxs'][3]])            
    bbb_obj.pick_nD(par[P['idxs'][3]:P['idxs'][4]])

    try: 
        bands, gal_Fnu=  gal_obj.get_fluxes(gal_obj.matched_parkeys)
        _, sb_Fnu= sb_obj.get_fluxes(sb_obj.matched_parkeys)
        _, bbb_Fnu = bbb_obj.get_fluxes(bbb_obj.matched_parkeys)
        _, tor_Fnu= tor_obj.get_fluxes(tor_obj.matched_parkeys)

    except ValueError:
         print ('Error: Dictionary does not contain some values')

    if models.settings['BBB'] !='R06':  #The other accretion disk models have a different normalization (only during the exploration of the parameters space)
        bbb_Fnu = bbb_Fnu/ (4*np.pi*dlum**2)
        BB=0

    # Total SED sum
    #--------------------------------------------------------------------
    lum = 10**(SB)* sb_Fnu  + 10**(BB)*bbb_Fnu + 10**(GA)*gal_Fnu +10**(TO) *tor_Fnu

    if models.settings['RADIO'] == True:  #Include the 5th component, if radio data is available
        lum += 10**(RAD)*agnrad_Fnu
    #--------------------------------------------------------------------    
    lum = lum.reshape((np.size(lum),))
    return lum, bands


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
        p0[:, i] = 0.5*(P['min'][i] + P['max'][i]) + (2* p0[:, i] - 1) *0.00001 
    
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
        
        p1[:, i] = p + p1[:, i]*0.00001 

    return p1




