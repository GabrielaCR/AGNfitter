
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
import pickle
import MODEL_AGNfitter as model
from DATA_AGNfitter import DATA


def Pdict (data):

    """Constructs a dictionary P with the name of all model parameters as keys. 
    The values are tuples with the same length, with the parameter limits. 

    ## inputs:
    - data object
    """
    P = adict()
    
    P.names ='tau', 'age', 'Nh', 'irlum' , 'SB', 'BB', 'GA', 'TO', 'EBVbbb', 'EBVgal'
    P.min = 0., 5., 21., 7., 0., 0., 0., 0., -0.1, -0.1
    P.max = 15, np.log10(model.maximal_age(data.z)), 25, 15, 10., 10., 10., 10., 1.0, 1.0

    Npar = len(P.names)

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
        if not (P.min[i] < p < P.max[i]):
            return -np.inf

    #2. Prior on the luminosity
    B_band_expected, B_band_thispoint = galaxy_Lumfct_prior( z, dlum, bands, gal_Fnu)    # Bband expectations
    #if Bband magnitude in this trial is brighter than expected by the luminosity function, dont accept this one
    if B_band_thispoint < (B_band_expected - 5):#2.5):
        return -np.inf

    return 0.



def ln_likelihood(x, y, ysigma, z, y_model):

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

    resid = [(y[i] - y_model[i])/ysigma[i] for i in x_valid]


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

    y_model, bands, gal_Fnu = ymodel(data.nus,data.z,data.dictkey_arrays,data.dict_modelfluxes,*pars)

    lnp = ln_prior(data.z, data.dlum, bands,gal_Fnu, P, pars)

    if np.isfinite(lnp):    
        posterior = lnp + ln_likelihood(data.nus,data.fluxes,data.fluxerrs, data.z, y_model)     
        return posterior
    return -np.inf



"""------------------------------------
CONSTRUCT TOTAL MODEL 
------------------------------------"""

def ymodel(data_nus, z, dictkey_arrays, dict_modelfluxes, *par):

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
    STARBURSTFdict , BBBFdict, GALAXYFdict, TORUSFdict,_,_,_,_,_= dict_modelfluxes
    gal_do, irlum_dict, nh_dict, BBebv_dict,_= dictkey_arrays

    # Call MCMC-parameter values 
    tau, agelog, nh, irlum, SB ,BB, GA,TO, BBebv, GAebv= par[0:10]
    age = 10**agelog

    # Pick dictionary key-values, nearest to the MCMC- parameter values
    irlum_dct = model.pick_STARBURST_template(irlum, irlum_dict)
    nh_dct = model.pick_TORUS_template(nh, nh_dict)
    ebvbbb_dct = model.pick_BBB_template(BBebv, BBebv_dict)
    gal_do.nearest_par2dict(tau, age, GAebv)
    tau_dct, age_dct, ebvg_dct=gal_do.t, gal_do.a,gal_do.e

    # Call fluxes from dictionary using keys-values
    try: 
        bands, gal_Fnu = GALAXYFdict[tau_dct, age_dct,ebvg_dct]     
        bands, sb_Fnu= STARBURSTFdict[irlum_dct] 
        bands, bbb_Fnu = BBBFdict[ebvbbb_dct] 
        bands, tor_Fnu= TORUSFdict[nh_dct]
 
    except ValueError:
        print 'Error: Dictionary does not contain some values'

    # Renormalize to have similar amplitudes. Keep these fixed!
    sb_Fnu_norm = sb_Fnu.squeeze()/1e20    
    bbb_Fnu_norm = bbb_Fnu.squeeze()/1e60
    gal_Fnu_norm = gal_Fnu.squeeze()/1e18
    tor_Fnu_norm = tor_Fnu.squeeze()/1e-40


    # Total SED sum
    #--------------------------------------------------------------------

    lum = 10**(SB)* sb_Fnu_norm + 10**(BB)*bbb_Fnu_norm    \
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

    Npar = len(P.names) 
    p0 = np.random.uniform(size=(nwalkers, Npar))

    for i in range(Npar):

        p0[:, i] = 0.5*(P.max[i] + P.min[i]) + (2* p0[:, i] - 1) * (1)
        
    p0[:,8]= 0.1 + (2* p0[:, 8] - 1) * (0.001)
    
    return p0


def get_best_position(filename, nwalkers, P):

    """Returns the best positions after burn-in phases.
    ## inputs:
    - filename (str), nwalkers(int), P (dict)
    ## output:
    - list of positions of length (P.keys)"""

    Npar = len(P.names) 
    #all saved vectors    
    f = open(filename, 'rb')
    samples = pickle.load(f)
    f.close()

    #index for the largest likelihood     
    i = samples['lnprob'].ravel().argmax()
    #the values for the parameters at this index
    P.ml = samples['chain'].reshape(-1, Npar)[i]

    p1 = np.random.normal(size=(nwalkers, Npar))

    for i in range(Npar):
        p = P.ml[i]
        
        p1[:, i] = p + 0.00001 * p1[:, i]

    return p1







class adict(dict): 

    """ A dictionary with attribute-style access. It maps attribute
    access to the real dictionary.
    This class has been obtained from the Barak package by Neil Chrighton)
    """

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    # the following two methods allow pickling
    def __getstate__(self):
        """Prepare a state of pickling."""
        return self.__dict__.items()

    def __setstate__(self, items):
        """ Unpickle. """
        for key, val in items:
            self.__dict__[key] = val

    def __setitem__(self, key, value):
        return super(adict, self).__setitem__(key, value)

    def __getitem__(self, name):
        return super(adict, self).__getitem__(name)

    def __delitem__(self, name):
        return super(adict, self).__delitem__(name)

    def __setattr__(self, key, value):
        if hasattr(self, key):
            # make sure existing methods are not overwritten by new
            # keys.
            return super(adict, self).__setattr__(key, value)
        else:
            return super(adict, self).__setitem__(key, value)

    __getattr__ = __getitem__

    def copy(self):
        """ Return a copy of the attribute dictionary.

        Does not perform a deep copy
        """
        return adict(self)

