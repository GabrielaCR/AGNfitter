import numpy as np
from math import exp,pi, sqrt
import matplotlib.pyplot as plt
import time
from astropy.table import Table
from astropy.io import fits, ascii
import scipy
import astropy.constants as const
import astropy.units as u
import itertools


def PRIORS(data, models, P, *pars):

    modelsettings= models.settings

    _ , BBBFdict, GALAXYFdict, _,_,_,_,_,GALAXY_SFRdict, GALAXYatt_dict, STARBURST_LIRdict, _ = models.dict_modelfluxes
    gal_obj,sb_obj,tor_obj, bbb_obj = models.dictkey_arrays

    if len(bbb_obj.par_names)==1:
        GA, SB, TO, BB= pars[-4:]
    else:
        GA, SB, TO = pars[-3:]
        BB = 0  

    all_priors=[]

    if modelsettings['PRIOR_energy_balance'] == True:  
        
        prior= prior_energy_balance(GALAXYatt_dict, gal_obj, GA, STARBURST_LIRdict,sb_obj,SB)
        all_priors.append(prior)

    ### Informative priors to be added

    if modelsettings['PRIOR_galaxy_only']==True:  
        """
        """
        prior= prior_low_AGNfraction(data, models, P, *pars)
        all_priors.append(prior)

    if modelsettings['PRIOR_AGNfraction']==True:  
        """
        """
        t1= time.time()

        prior1= prior_AGNfraction(data, GALAXYFdict, gal_obj, GA, BBBFdict, bbb_obj, BB)
        prior2= prior_stellar_mass(GA)
        all_priors.append(prior1+prior2)
        #print('INTERNAL',prior, time.time()-t1)


    if modelsettings['XRAYS']==True:  

        prior=  prior_xrays(data, models, P, *pars)
        all_priors.append(a)

    if modelsettings['RADIO']==True:  
        a=0
        all_priors.append(a)

    final_prior= np.sum(np.array(all_priors))


    return final_prior


def prior_energy_balance(GALAXYatt_dict, gal_obj, GA, STARBURST_LIRdict,sb_obj,SB):

    Lgal_att = GALAXYatt_dict[gal_obj.matched_parkeys] * 10**(GA)
    Lsb_emit = STARBURST_LIRdict[sb_obj.matched_parkeys] * 10**(SB)

    if Lsb_emit < Lgal_att:
        return -np.inf
    else:
        return 0


def prior_AGNfraction(data, GALAXYFdict, gal_obj,GA, BBBFdict, bbb_obj, BB): 

    bands, gal_Fnu= GALAXYFdict[gal_obj.matched_parkeys]
    bands, bbb_Fnu = BBBFdict[bbb_obj.matched_parkeys] 

    gal_flux= gal_Fnu* 10**(GA)
    bbb_flux= bbb_Fnu* 10**(BB)

    """calculate 1500 Angstrom magnitude in the data and model"""

    data_flux_1500Angs = data.fluxes[(14.3 < (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 15.3 )]
    data_flux_1500Angs=data_flux_1500Angs[data_flux_1500Angs>0][-1]

    gal_flux_1500Angs = gal_flux[(14.3 < bands+np.log10(1+data.z)) & (bands+np.log10(1+data.z) < 15.3)][-1]
    bbb_flux_1500Angs = bbb_flux[(14.3 < bands+np.log10(1+data.z)) & (bands+np.log10(1+data.z) < 15.3)][-1]
  
    data_lum_1500Angs = data.lumfactor*data_flux_1500Angs
    abs_mag_data = 51.6 - 2.5 *np.log10(data_lum_1500Angs)

    characteristic_mag = -35.4 * (1+data.z)**0.524/(1+(1+data.z)**0.678) ### Expected UV magnitude from Parsa, Dunlop et al. 2014.
                                                                         ### these calculations are based on Hubble Ultra Deep Field (HUDF), CANDELS/GOODS-South,
                                                                         ### and UltraVISTA/COSMOS surveys data from z~ 2-4, and literature at lower redshifts.

    """define prior on agnfraction"""
    if BB ==0:
        bbb_flux_1500Angs = bbb_flux_1500Angs/(4*pi*(data.dlum)**2)   ##BB normalization
    AGNfrac1500 = np.log10(bbb_flux_1500Angs/gal_flux_1500Angs) 

    if abs_mag_data > (characteristic_mag-1.): ## if blue fluxes are fainter than 10 times the characteristic flux.
                                               ## asume galaxy dominates, unless data strongly prefers so.
        mu = -2.
        sigma = 2.
        prior_AGNfrac = Gaussian_prior(mu, sigma, AGNfrac1500)

    else:                                      ## if blue fluxes are equal or brighter than 10 times the characteristic flux.
                                               ## asume BBB is at least equal to galaxy or dominates.
        if AGNfrac1500<0:
            prior_AGNfrac=-np.inf
        else:
            ### Adding the prior knowledge on AGN fraction only valid for QSO.
            mu = 2
            sigma = 2.
            prior_AGNfrac = Gaussian_prior(mu, sigma, AGNfrac1500)

    return prior_AGNfrac


def prior_stellar_mass(GA):
    ### Adding the prior knowledge on stellar masses of host galaxies
    ### GA<3 corresponds to M* < 1e9 Msun/yr
    mu_GA =4.5
    sigma_GA =1.5
    prior_GA = Gaussian_prior(mu_GA, sigma_GA, GA)

    return prior_GA 


# def prior_xrays(data, models, P, *pars):

# 	def alpha_OX(log_L2kev):
# 		""" Relation between accretion disk intrinsic luminosity at 2500 Angstrom 
# 			and X-rays at 2 keV .
# 			Lusso&Risaliti +16 gives beta=[0.6-0.65], gamma=[7-8]"""
# 		beta= 0.6 
# 		gamma= 7.5 
# 		log_L2500A_alphaox= (log_L2kev-beta)/gamma

# 		return log_L2500A_alphaox

#     _ , BBBFdict, _, _,_,_,_,_, _, _, _, _ = models.dict_modelfluxes

#     if len(bbb_obj.par_names)==1:
#         GA, SB, TO, BB= pars[-4:]
#     else:
#         GA, SB, TO = pars[-3:]

# 	all_bbb_nus, bbb_Fnus_dered = BBBFdict['0.0']
#     bbb_flux_dered= bbb_Fnus_dered* 10**(BB)

#     """Calculate 2kev (10**17.684 Hz) and 2500 Angstrom (10**15.06 Hz) magnitude in the data and model"""

#     flux_2kev = data.fluxes[( 17.60 < (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 17.75 )]

#     bbb_flux_dered_2500Angs = data.fluxes[(15.04< (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 15.15 )]
#     if len(bbb_flux_dered_2500Angs)>1:
#         bbb_flux_dered_2500Angs = bbb_flux_dered_2500Angs[0]
#     lumfactor = (4. * pi * data.dlum**2.)
#     log_L2500A_data_dered = np.log10(lumfactor*bbb_flux_dered_2500Angs)

#     log_L2500A_model= alpha_OX(log_L2kev)

#     ratio_alpha0x_data= log_L2500A_data_dered - log_L2500A_model

#     """Define prior"""
#     mu= 0
#     sigma= 0.4
#     prior_Xrays= Gaussian_prior(mu, sigma, ratio_alpha0x_data)

#     return prior_Xrays


def prior_low_AGNfraction(data, models, P, *pars):

    _ , BBBFdict, GALAXYFdict, _,_,_,_,_, _, GALAXYatt_dict, _, _ = models.dict_modelfluxes
    gal_obj,_,_, bbb_obj = models.dictkey_arrays

    if len(bbb_obj.par_names)==1:
        GA, SB, TO, BB= pars[-4:]
    else:
        GA, SB, TO = pars[-3:]

    bands, gal_Fnu= GALAXYFdict[gal_obj.matched_parkeys]
    bands, bbb_Fnu = BBBFdict[bbb_obj.matched_parkeys] 

    gal_flux= gal_Fnu* 10**(GA)
    bbb_flux= bbb_Fnu* 10**(BB)

    """calculate 1500 Angstrom magnitude in the data and model"""

    data_flux_1500Angs = data.fluxes[(15. < (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 15.3 )]
    gal_flux_1500Angs = gal_flux[(14.7 < bands+np.log10(1+data.z)) & (bands+np.log10(1+data.z) < 15.3)]
    bbb_flux_1500Angs = bbb_flux[(14.7 < bands+np.log10(1+data.z)) & (bands+np.log10(1+data.z) < 15.3)]

    if data_flux_1500Angs[-1]>0:### Check if it's not a non-detection (-99)
        data_flux_1500Angs=data_flux_1500Angs
    else:
        data_flux_trial= data.fluxes[(14. < (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 15.3 )]
        data_flux_1500Angs=data_flux_trial[data_flux_trial>0]

    if len(data_flux_1500Angs)>1:
        gal_flux_1500Angs = gal_flux_1500Angs[-1]
        bbb_flux_1500Angs = bbb_flux_1500Angs[-1]
        data_flux_1500Angs = data_flux_1500Angs[data_flux_1500Angs>0][-1]

    lumfactor = (4. * pi * data.dlum**2.)
    data_lum_1500Angs = lumfactor*data_flux_1500Angs
    abs_mag_data = 51.6 - 2.5 *np.log10(data_lum_1500Angs)

    """ Expected UV magnitude from Parsa, Dunlop et al. 2016.
    These calculations are based on Hubble Ultra Deep Field (HUDF), CANDELS/GOODS-South,
    and UltraVISTA/COSMOS surveys data from z~ 2-4, and literature at lower redshifts."""
    characteristic_mag = -35.4 * (1+data.z)**0.524/(1+(1+data.z)**0.678)

    """Setting-up prior"""
    AGNfrac1500 = np.log10(bbb_flux_1500Angs/gal_flux_1500Angs)
    if len(AGNfrac1500)>1:
        AGNfrac1500=AGNfrac1500[0]
    if abs_mag_data > (characteristic_mag-3.): # If UV luminosity  is below the characteristic galaxy luminosity at that given redshifts
                                          # the luminosity is preferable fitted by the stellar component rather than the AGN,
                                          # unless the data strongly prefers it.
        mu = -2.
        sigma = 0.5
        prior_AGNfrac = Gaussian_prior(mu, sigma, AGNfrac1500)
        #prior_AGNfrac = Clipped_Gaussian_prior(mu, sigma, -5, 0.2, AGNfrac1500)        
        #prior_AGNfrac = Clipped_TophatAndGaussian_prior(mu, sigma, -5, 0.2, AGNfrac1500)
    else: ##type2
        mu = -2.
        sigma = 2
        prior_AGNfrac = Gaussian_prior(mu, sigma, AGNfrac1500)

    return prior_AGNfrac


def Gaussian_prior(mu, sigma, par):

	return  np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5 * (par - mu)**2/sigma**2  	

def Clipped_Gaussian_prior(mu, sigma, clipmin, clipmax, par):

    if par>clipmin and par<clipmax:
        return  np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5 * (par - mu)**2/sigma**2  
    else:
        return -np.inf

def Clipped_TophatAndGaussian_prior(mu, sigma, clipmin, clipmax, par):

    if par>=mu and par<clipmax:
        return  np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5 * (par - mu)**2/sigma**2  
    elif par<mu and par>clipmin :
        return 0
    else:
    	return -np.inf

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


def maximal_age(z):
    """
    Maximal possible age for galaxy, set by the age of the Universe at the time of observation.
    """

    z = np.double(z)
    #Cosmological Constants    
    O_m = 0.266
    O_r =  0.
    O_k= 0.
    O_L = 1. - O_m
    H_0 = 74.3 #km/s/Mpc
    H_sec = H_0 / 3.0857e19 
    secondsinyear = 31556926
    ageoftheuniverse = 13.798e9

    # Equation for the time elapsed since z and now

    a = 1/(1+z)
    E = O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L
    integrand = lambda z : 1 / (1+z)     / sqrt(  O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L  )        

    #Integration
    z_obs = z
    z_cmb = 1089  
    z_now = 0


    integral, error = scipy.integrate.quad( integrand , z_obs, z_cmb) #
    
    #t = ageoftheuniverse - (integral * (1 / H_sec) / secondsinyear)
    t = (integral * (1 / H_sec) / secondsinyear) - 1e9

    return t
