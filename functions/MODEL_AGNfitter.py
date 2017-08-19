

"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         MODEL_AGNfitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains all functions which are needed to construct the total model of AGN. 
The functions here translate the parameter space points into total fluxes dependin on the models chosen.

Functions contained here are the following:

STARBURST_nf
BBB_nf
GALAXY_nf
TORUS_nf

"""

import numpy as np
from math import exp,pi, sqrt
import matplotlib.pyplot as plt
import time
import cPickle
from astropy.table import Table
from astropy.io import fits
from scipy.integrate  import quad, trapz
import astropy.constants as const
import astropy.units as u
import itertools
                     
"""
ADDING NEW MODELS
------------------
- Add the model similar to the examples given below, 
    if modelsettings['COMPONENT']=='MYMODEL':
        ....
        my_parameternames= ['par1', 'par2']
        MYMODELdict['par1v','par2v'] = log10(nus[Hz]), Lnu 
        ....
        return MYMODELdict_4plot,  my_parameternames
- Write the name you gave 'MYMODEL' to your settings file for the respectively. 

- Other changes to be done manually (unfortuantely):
  (Due to computational time reasons not possible yet to do this automatically)

    (1) Go to fct ymodel in PARAMETERSPACE_AGNfitter.py and set (eg. for the galaxy model)
        gal_obj.pick_1D (if your model has one parameter)
        or 
        gal_obj.pick_nD (if your model has more parameters).

    (2) In the same function if your model has more than one parameter, change:
        GALAXYFdict[gal_obj.matched_parkeys] (model of one parameter)
        or
        GALAXYFdict[tuple(gal_obj.matched_parkeys)] (model of more than one parameters

    (3) and (4) Go to PLOTandWRITE_AGNfitter.py, method fluxes of class FLUXES_ARRAYS, 
    and do the same two changes as above.


"""
                                                         

def GALAXY(path, modelsettings):

    if modelsettings['GALAXY']=='BC03':


        GALAXYFdict_4plot = dict()
        GALAXY_SFRdict = dict()
        ## Call object containing all galaxy models     
        BC03dict = cPickle.load(file(path + 'models/GALAXY/BC03_275seds.pickle', 'rb'))    

        ## specify the sizes of the array of parameter values: Here two parameters
        tau_array = BC03dict['tau-values']
        age_array = BC03dict['age-values']
        ebvgal_array = np.array(np.arange(0.,100.,5.)/100)

        ## produce all combinations of parameter values (indices)
        _, ageidx, tauidx, _, _,_ =  np.shape(BC03dict['SED'])
        idxs = [np.arange(ageidx), np.arange(tauidx), np.arange(len(ebvgal_array))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                agei=c[0]
                taui=c[1]
                ebvi=c[2]
                gal_wl, gal_Fwl =  BC03dict['wavelength'],BC03dict['SED'][:,agei,taui,:,:,:].squeeze()
                gal_nus= gal_wl.to(u.Hz, equivalencies=u.spectral())[::-1]#invert
                gal_Fnu= (gal_Fwl * 3.34e-19 * gal_wl**2.)[::-1]  
                gal_SFR= BC03dict['SFR'][:,agei,taui,:,:].squeeze()
                GALAXY_SFRdict[str(tau_array.value[taui]),str(age_array.value[agei])] = gal_SFR         
                gal_nu, gal_Fnu_red = GALAXYred_Calzetti(gal_nus.value[0:len(gal_nus):3], gal_Fnu.value[0:len(gal_nus):3], ebvgal_array[ebvi])                    
                GALAXYFdict_4plot[str(tau_array.value[taui]),str(age_array.value[agei]), str(ebvgal_array[ebvi])] = \
                                                                                        np.log10(gal_nu), gal_Fnu_red        

        ## Name the parameters that compose the keys of the dictionary: GALAXYFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['tau', 'age','EBVgal']

        return  GALAXYFdict_4plot, GALAXY_SFRdict, parameters_names


def STARBURST(path, modelsettings):

    if modelsettings['STARBURST']=='DH02_CE01':

        STARBURSTFdict_4plot = dict()

        #Call object containing all starburst models     
        DH02CE01dict = cPickle.load(file(path + 'models/STARBURST/DH02_CE01.pickle', 'rb')) 
        irlumidx = len(DH02CE01dict['SED'])

        #Construct dictionaries 
        for irlumi in range(irlumidx):
            sb_nu0, sb_Fnu0 = DH02CE01dict['wavelength'][irlumi], DH02CE01dict['SED'][irlumi].squeeze()
            print sb_nu0
            STARBURSTFdict_4plot[str(DH02CE01dict['irlum-values'][irlumi])] = sb_nu0, sb_Fnu0

        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['irlum']

        return STARBURSTFdict_4plot, parameters_names

    elif modelsettings['STARBURST']=='S17':

        STARBURSTFdict_4plot = dict()

        #Call object containing all starburst models     
        dusttable = Table.read(path + 'models/STARBURST/s17_dust.fits')
        pahstable = Table.read(path + 'models/STARBURST/s17_pah.fits')
        
        Dwl, DnuLnu = dusttable['LAM'],dusttable['SED'] #micron, Lsun
        Pwl, PnuLnu = pahstable['LAM'],pahstable['SED'] #micron, Lsun
        Tdust = np.array(dusttable['TDUST'])[0] #K
        fracPAH = np.arange(1.5,5.0,0.5)/100
        idxs=[np.arange(len(Tdust)), np.arange(len(fracPAH))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        Dnu= (Dwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        Pnu= (Pwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        DLnu= np.array(DnuLnu[0])/Dnu *1e-6 #* u.Lsun.to(u.W)
        PLnu=np.array(PnuLnu[0])/Pnu *1e-6#* u.Lsun.to(u.W)


        #Construct dictionaries 
        for c in par_idxs_combinations:
            t=c[0]
            fp=c[1]

            sb_nu0 = np.array(Dnu[t,:])[::-1]
            sb_Fnu0 = np.array( (1-fracPAH[fp]) * DLnu[t,:] + (fracPAH[fp]) * PLnu[t,:])[::-1]
            # print type(Tdust)#[t]
            # print (Dnu[t,:])
            # print np.shape(sb_Fnu0), sb_Fnu0

            STARBURSTFdict_4plot[str(Tdust[t]), str(fracPAH[fp])] = np.log10(sb_nu0), sb_Fnu0

        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['Tdust', 'fracPAH']

        return STARBURSTFdict_4plot, parameters_names


def BBB(path, modelsettings):

    if modelsettings['BBB']=='R06':

        BBBFdict_4plot = dict()
        R06dict = cPickle.load(file(path + 'models/BBB/R06.pickle', 'rb')) 
        #R06dict = cPickle.load(file(path + 'models/BBB/richards.pickle', 'rb')) 
        parameters_names =['EBVbbb']
        ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)

        bbb_nu, bbb_Fnu = R06dict['wavelength'], R06dict['SED'].squeeze()
        #bbb_nu, bbb_Fnu = R06dict.wave, R06dict.SED.squeeze()
        
        #Construct dictionaries
        for EBV_bbb in ebvbbb_array:
            bbb_nu0, bbb_Fnu_red = BBBred_Prevot(bbb_nu, bbb_Fnu, EBV_bbb )
            BBBFdict_4plot[str(EBV_bbb)] =bbb_nu0, bbb_Fnu_red

        ## Name the parameters that compose the keys of the dictionary: BBFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['EBVbbb']

        return BBBFdict_4plot, parameters_names


def TORUS(path, modelsettings):

    if modelsettings['TORUS']=='S04':    

        TORUSFdict_4plot  = dict()

        #Call object containing all torus models     
        S04dict = cPickle.load(file(path + 'models/TORUS/S04.pickle', 'rb')) 
        parameters_names = ['Nh']
        nhidx=len(S04dict['SED'])
        #Construct dictionaries 
        for nhi in range(nhidx):

            tor_nu0, tor_Fnu0 = S04dict['wavelength'][nhi], S04dict['SED'][nhi].squeeze()
            TORUSFdict_4plot[str(S04dict['Nh-values'][nhi])] = tor_nu0, tor_Fnu0

        ## Name the parameters that compose the keys of the dictionary: TORUSFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.   
        parameters_names = ['Nh']

        return TORUSFdict_4plot, parameters_names



#==============================
# MAXIMAL POSSIBLE AGE FOR GALAXY MODEL
#==============================


def maximal_age(z):

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
    z_cmb = 1089 #As Beta (not cmb). But 1089 (cmb) would be the exagerated maximun possible redshift for the birth 
    z_now = 0


    integral, error = quad( integrand , z_obs, z_cmb) #
    
    #t = ageoftheuniverse - (integral * (1 / H_sec) / secondsinyear)
    t = (integral * (1 / H_sec)) / secondsinyear

    return t



"""===================================================
Reddening functions    
==================================================="""


def BBBred_Prevot(bbb_x, bbb_y, BBebv ):

    """
    
    ## input:

    ## output:

    """
    #Application of reddening - reading E(B-V) from MCMC sampler
    RV= 2.72

    #converting freq to wavelenght, to be able to use prevots function instead on simple linera interpolation 
    redd_x =  2.998 * 1e10 / (10**(bbb_x)* 1e-8)
    redd_x= redd_x[::-1]

    #    Define prevots function for the reddenin law redd_k    
    def function_prevot(x, RV):
           y=1.39*pow((pow(10.,-4.)*x),-1.2)-0.38 ;
           return y 

    bbb_k = function_prevot(redd_x, RV)

    bbb_k= bbb_k[::-1]

    bbb_Lnu_red = bbb_y * 10**(-0.4 * bbb_k * BBebv)

    return bbb_x, bbb_Lnu_red


def GALAXYred_Calzetti(gal_nu, gal_Fnu,GAebv):

    """
    This function computes the effect of reddening in the galaxy template (Calzetti law)

    ## input:
    -frequencies in log nu
    - Fluxes in Fnu
    - the reddening value E(B-V)_gal
    ## output:

    """
    RV = 4.05        

    c =2.998 * 1e8 
    gal_lambda_m = c / gal_nu * 1e6#in um 
    wl = gal_lambda_m[::-1]  #invert for lambda
    k = np.zeros(len(wl))

    w0 = [wl <= 0.12]
    w1 = [wl < 0.63]
    w2 = [wl >= 0.63]

    x1 = np.argmin(np.abs(wl - 0.12))
    x2 = np.argmin(np.abs(wl - 0.125))

    k[w2] = 2.659 * (-1.857 + 1.040 /wl[w2])+RV
    k[w1] = 2.659 * (-2.156 + (1.509/wl[w1]) - (0.198/wl[w1]**2) + (0.011/wl[w1]**3))+RV
    k[w0] = k[x1] + ((wl[w0] - 0.12) * (k[x1] - k[x2]) / (wl[x1] - wl[x2])) +RV


    gal_k= k[::-1] #invert for nus
    gal_Fnu_red = gal_Fnu* 10**(-0.4 * gal_k * GAebv)

    return gal_nu, gal_Fnu_red


Angstrom = 1e10

def z2Dlum(z):

    """
    Calculate luminosity distance from redshift.
    """
    
    #Cosmo Constants
    
    O_m = 0.266
    O_r =  0.
    O_k= 0.
    O_L = 1. - O_m
    H_0 = 70. #km/s/Mpc
    H_sec = H_0 / 3.0857e19 
    # equation

    a = 1/(1+z)
    E = O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L
    integrand = lambda z : 1 / sqrt(O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L)    

    #integration

    z_obs = z
    z_now = 0

    c_cm = 2.997e10

    
    integral = quad( integrand , z_now, z_obs)    
    dlum_cm = (1+z)*c_cm/ H_sec * integral[0] 
    dlum_Mpc = dlum_cm/3.08567758e24

    return dlum_cm
   

"""---------------------------------------------
             COMPUTED QUANTITIES
-----------------------------------------------"""

def stellar_info(chain, data):

    """
    computes stellar masses and SFRs
    """

    gal_obj,_,_,_ = data.dictkey_arrays
    _,_,_,_,_,_,_,_,SFRdict,_= data.dict_modelfluxes
    #relevanta parameters form the MCMC chain
    tau_mcmc = chain[:,0]     
    age_mcmc = chain[:,1] 
    GA = chain[:, -4] - 18. ## 1e18 is the common normalization factor used in parspace.ymodel 
                            ## in order to have comparable NORMfactors    

    z = data.z
    distance = z2Dlum(z)

    #constants
    solarlum = const.L_sun.to(u.erg/u.second) #3.839e33
    solarmass = const.M_sun

    Mstar_list=[]
    SFR_list=[]

    print len(tau_mcmc), len(GA)
    for i in range (len (tau_mcmc)):        
        N = 10**GA[i]* 4* pi* distance**2 / (solarlum.value)/ (1+z)

        gal_obj.pick_nD(tuple([tau_mcmc[i], age_mcmc[i], 0.]))
        tau_dct, age_dct, ebvg_dct=gal_obj.matched_parkeys
        SFR_mcmc =SFRdict[tau_dct, age_dct]

        # Calculate Mstar. BC03 templates are normalized to M* = 1 M_sun. 
        # Thanks to Kenneth Duncan, and his python version of BC03, smpy
        Mstar = np.log10(N * 1) 
        #Calculate SFR. output is in [Msun/yr]. 
        SFR = N * SFR_mcmc
        SFR_list.append(SFR.value)    
        Mstar_list.append(Mstar)    

    return np.array(Mstar_list) , np.array(SFR_list)


def stellar_info_array(chain_flat, data, Nthin_compute):

    """
    computes arrays of stellar masses and SFRs
    """

    Ns, Npar = np.shape(chain_flat)  
    chain_thinned = chain_flat[0:Ns:int(Ns/Nthin_compute),:]

    Mstar, SFR = stellar_info(chain_thinned, data)
    Mstar_list = []
    SFR_list = []

    for i in range(Nthin_compute):
        for j in range(int(Ns/Nthin_compute)):
            Mstar_list.append(Mstar[i])
            SFR_list.append(SFR[i])

    Mstar1 = np.array(Mstar_list)    
    SFR1 = np.array(SFR_list)
    return Mstar1, SFR1



def sfr_IR(logL_IR):
    #calculate SFR in solar M per year 

    #for an array ofluminosities
    if len(logL_IR)>1:
        SFR_IR_list =[]

        for i in range(len(logL_IR)):
            SFR = 3.88e-44* (10**logL_IR[i])
            SFR_IR_list.append(SFR)
        SFR_IR_array = np.array(SFR_IR_list)
        return SFR_IR_array
    #or just for one luminosity
    else:        
        SFR = 3.88e-44* (10**logL_IR)
        return SFR


