import numpy as np
from math import pi, sqrt
import time
import scipy


def PRIORS(data, models, P, *pars):

    modelsettings= models.settings
    MD = models.dict_modelfluxes
    gal_obj,sb_obj,tor_obj, bbb_obj, agnrad_obj = models.dictkey_arrays

    if modelsettings['BBB']=='R06':
        if models.settings['RADIO'] == True:
            GA, SB, TO, BB, RAD = pars[-5:]
        else:
            GA, SB, TO, BB= pars[-4:]
    else:
        BB = 0                      #If accretion disk model is different from R06, the normalization is different
        if models.settings['RADIO'] == True:
            GA, SB, TO, RAD = pars[-4:]
        else:
            GA, SB, TO = pars[-3:]

    all_priors=[]

    if (modelsettings['PRIOR_energy_balance'] == 'Flexible') or (modelsettings['PRIOR_energy_balance'] == 'Restrictive'):  
        """
        This prior promotes starburst emission consistent with galaxy attenuated emission. The flexible prior only impose a lower limit
        for the luminosity of the cold dust, while the restrictive promotes models in which both emissioons are the same.
        """
        prior= prior_energy_balance(data, MD.GALAXYatt_dict, MD.GALAXYFdict, gal_obj, GA, MD.STARBURST_LIRdict,sb_obj,SB, models)
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
        prior1= prior_AGNfraction(data, MD.GALAXYFdict, gal_obj, GA, MD.BBBFdict, bbb_obj, BB)
        prior2= prior_stellar_mass(GA)
        all_priors.append(prior1 + prior2) 

    if modelsettings['PRIOR_midIR_UV']==True:  
        """
        """
        prior_IR_UV= prior_midIR_UV(data, MD.BBBFdict, bbb_obj, BB, MD.TORUSFdict, tor_obj, TO, models)
        all_priors.append(prior_IR_UV) 


    if modelsettings['RADIO']==True:  
        # This prior gives predominance to Synchrotron more than Starburst emission in IR if the IR data can be explained by a simple power law 
        # extended from radio data available
        prior_radio = prior_IR_SYNfraction(data, MD.STARBURSTFdict, sb_obj, SB, MD.AGN_RADFdict, agnrad_obj, RAD, models)
        all_priors.append(prior_radio)

    if modelsettings['XRAYS']== 'Prior_UV': 
        # This prior promotes accretion disk models consistent with Xrays data and the alpha_ox correlation by Just et al. 2007
        prior_L2500_alphaox = prior_UV_xrays(data, MD.BBBFdict, bbb_obj, BB, models)
        all_priors.append(prior_L2500_alphaox)

    if modelsettings['XRAYS']== 'Prior_midIR': 
        # This prior promotes torus models consistent with Xrays data and the mid-IR-Xray correlation by Stern 2015
        prior_L6microns = prior_IR_XRays(data, MD.TORUSFdict, tor_obj, TO, models)
        all_priors.append(prior_L6microns)

    final_prior= np.sum(np.array(all_priors))

    return final_prior


def prior_energy_balance(data, GALAXYatt_dict, GALAXYFdict, gal_obj, GA, STARBURST_LIRdict,sb_obj,SB, models):

    if gal_obj.par_types[-1] == 'grid':
        Lgal_att = GALAXYatt_dict[tuple(gal_obj.matched_parkeys_grid)] * 10**(GA)

    elif gal_obj.par_types[-1] == 'free':
        bands, gal_Fnu= GALAXYFdict[tuple(gal_obj.matched_parkeys_grid)]        #frequencies in log
        fcts=gal_obj.functions()
        f=fcts[gal_obj.functionidxs[0]]
        rest_bands = bands + np.log10((1+data.z))                               #Pass to rest frame
        bandsf, Fnuf = f(10**rest_bands, gal_Fnu*1e18, gal_obj.matched_parkeys[-1])  #bandsf not in log form, apply reddening
        gal_nu, gal_Fnu_red = bandsf/(1+data.z), Fnuf                           #Pass to observed frame
        gal_Fnu_int = scipy.integrate.trapz(gal_Fnu*3.826e33, x=gal_nu)          
        gal_Fnured_int = scipy.integrate.trapz(gal_Fnu_red*3.826e33, x=gal_nu)
        gal_att_int = gal_Fnu_int - gal_Fnured_int
        Lgal_att = abs(gal_att_int * 10**(GA))                                       #Calculate the attenuated luminosity

    Lsb_emit = STARBURST_LIRdict[sb_obj.matched_parkeys] * 10**(SB)

    if Lsb_emit < Lgal_att:
        return -np.inf
    elif (Lsb_emit >= Lgal_att) and (models.settings['PRIOR_energy_balance'] == 'Flexible'):
        return 0
    elif (Lsb_emit >= Lgal_att) and (models.settings['PRIOR_energy_balance'] == 'Restrictive'):
        frac_SB_attGal = np.log10(Lsb_emit/Lgal_att)
        mu = 0
        sigma = 0.1
        prior_frac = Gaussian_prior(mu, sigma, frac_SB_attGal)
        return prior_frac


def prior_AGNfraction(data, GALAXYFdict, gal_obj,GA, BBBFdict, bbb_obj, BB): 

    if gal_obj.par_types[-1] == 'grid':
        bands, gal_Fnu= GALAXYFdict[tuple(gal_obj.matched_parkeys_grid)]
    elif gal_obj.par_types[-1] == 'free':
        bands, gal_Fnu= GALAXYFdict[tuple(gal_obj.matched_parkeys_grid)]
        fcts=gal_obj.functions()
        f=fcts[gal_obj.functionidxs[0]]
        rest_bands = bands + np.log10((1+data.z))                               #Pass to rest frame
        bandsf, Fnuf = f(10**rest_bands, gal_Fnu, gal_obj.matched_parkeys[-1])  #Apply galaxy reddening
        bandsf = np.log10(bandsf) - np.log10((1+data.z))                        #Come back to observed frame
        bands, gal_Fnu = bandsf, Fnuf

    if bbb_obj.par_types[-2: ] == ['free', 'free'] and bbb_obj.par_names[-2: ] == ['EBVbbb', 'alphaScat']: 
        fcts=bbb_obj.functions()
        f=fcts[bbb_obj.functionidxs[0]]
        bands, Fnu = bbb_obj.modelsdict[tuple(bbb_obj.matched_parkeys_grid)]
        rest_bands = bands + np.log10((1+data.z))                               #Rest frame frequency
        bandsf0, Fnuf0 = f(rest_bands[rest_bands < 16.685], Fnu[rest_bands < 16.685], bbb_obj.matched_parkeys[-2])  #Apply reddening
        bandsf =  np.concatenate((bandsf0, rest_bands[rest_bands >= 16.685])) - np.log10((1+data.z))            #Come back to observed frame 
        Fnuf = np.concatenate((Fnuf0, Fnu[rest_bands >= 16.685]*10**bbb_obj.matched_parkeys[-1]))               #Add the effect of scatter in UV-Xray correlation
        bands, bbb_Fnu = bandsf, Fnuf

    elif bbb_obj.par_types[-3: ] == ['free', 'free', 'grid'] and bbb_obj.par_names[-3: ] == ['EBVbbb', 'alphaScat', 'Gamma']: 
        fcts=bbb_obj.functions()
        f=fcts[bbb_obj.functionidxs[0]]
        bands, Fnu = bbb_obj.modelsdict[tuple(bbb_obj.matched_parkeys_grid)]
        rest_bands = bands + np.log10((1+data.z))                               #Rest frame frequency
        bandsf0, Fnuf0 = f(rest_bands[rest_bands < 16.685], Fnu[rest_bands < 16.685], bbb_obj.matched_parkeys[-3])  #Apply reddening
        bandsf =  np.concatenate((bandsf0, rest_bands[rest_bands >= 16.685])) - np.log10((1+data.z))            #Come back to observed frame 
        Fnuf = np.concatenate((Fnuf0, Fnu[rest_bands >= 16.685]*10**bbb_obj.matched_parkeys[-2]))               #Add the effect of scatter in UV-Xray correlation
        bands, bbb_Fnu = bandsf, Fnuf

    elif bbb_obj.par_types[-1] == 'free' and bbb_obj.par_names[-1] == 'EBVbbb':  #This is for the case of EBVbbb == free and no x-rays
        fcts=bbb_obj.functions()
        f=fcts[bbb_obj.functionidxs[0]]
        if type(bbb_obj.matched_parkeys_grid) != list:                          #If model is R06, there isn't a list of parameters so tuple don't work well
            bands, Fnu = bbb_obj.modelsdict[bbb_obj.matched_parkeys_grid] 
            bbb_obj.matched_parkeys = [bbb_obj.matched_parkeys]
        else:
            bands, Fnu = bbb_obj.modelsdict[tuple(bbb_obj.matched_parkeys_grid)] 
        rest_bands = bands + np.log10((1+data.z))                               #Rest frame frequency
        bandsf, Fnuf = f(rest_bands, Fnu, bbb_obj.matched_parkeys[-1])          #Apply reddening
        bandsf = bandsf - np.log10((1+data.z))                                  #Come back to observed frame 
        bands, bbb_Fnu = bandsf, Fnuf  
    else:                                                                       #When all parameters are grid
        bands, bbb_Fnu = BBBFdict[bbb_obj.matched_parkeys] 


    gal_flux= gal_Fnu* 10**(GA)
    bbb_flux= bbb_Fnu* 10**(BB)

    """calculate 1500 Angstrom magnitude in the data and model"""

    data_flux_1500Angs = data.fluxes[(14.3 < (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 15.3 )]
    data_flux_1500Angs = data_flux_1500Angs[data_flux_1500Angs>0][-1]

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



def prior_IR_SYNfraction(data, STARBURSTFdict, sb_obj, SB, AGN_RADFdict, agnrad_obj, RAD, models): 

    bands, sb_Fnu = STARBURSTFdict[sb_obj.matched_parkeys] 
    if (agnrad_obj.pars_modelkeys != ['-99.9']).all() :   #If the model have fitting parameters
        bands, syn_Fnu = AGN_RADFdict[agnrad_obj.matched_parkeys]
    else:
        bands, syn_Fnu = AGN_RADFdict['-99.9']          #If the model have fix parameters
    sb_flux= sb_Fnu* 10**(SB)
    syn_flux = syn_Fnu* 10**(RAD)

    """calculate IR flux from synchrotron and data"""

    data_flux_rad = data.fluxes[((8-np.log10(1+data.z)) < data.nus) & (data.nus < (10-np.log10(1+data.z)) )]   #0.1-10 GHz rest frame --> to observed frame
    if True not in (data_flux_rad>0):                                                                          # Not valid data == not prior information
        return 0
    data_flux_rad = data_flux_rad[data_flux_rad>0][-1]                                                         #We choose the highest frequency data
    data_nu_rad = data.nus[data.fluxes == data_flux_rad][0]                                                    #its frequency


    data_flux_IR = data.fluxes[((12.2-np.log10(1+data.z)) < data.nus) & (data.nus < (12.8-np.log10(1+data.z)))] #peak cold dust spectrum rest frame
    if True not in (data_flux_IR>0):
        return 0
    data_flux_IR = max(data_flux_IR[data_flux_IR>0])                                                            #We choose the highest flux
    data_nu_IR = data.nus[data.fluxes == data_flux_IR][0]                                                       #its frequency

    syn_exp_IR = data_flux_rad*((10**data_nu_IR/10**data_nu_rad)**(-0.75))        #Expected IR flux if synchrotron emission dominates (single power law)

    """define prior on SYNfrac_IR"""
    sb_flux_IR = sb_flux[np.argmin(np.abs(bands-data_nu_IR))]
    syn_flux_IR = syn_flux[np.argmin(np.abs(bands-data_nu_IR))]  
    SYNfrac_IR = np.log10(syn_flux_IR/sb_flux_IR)

    if (data_flux_IR/syn_exp_IR) < 2 :                    ## IR flux is near to the value estimated from the synchrotron power law, asume RAD dominates      
        mu = 2.
        sigma = 2.                                                           
        prior_SYNfrac = Gaussian_prior(mu, sigma, SYNfrac_IR)

    elif (data_flux_IR/syn_exp_IR) > 2 :                 ## if IR flux from synchrotron power law is very low, STARBURST dominates
        mu = -2
        sigma = 2.
        prior_SYNfrac = Gaussian_prior(mu, sigma, SYNfrac_IR)                                 
 
    return prior_SYNfrac

def prior_UV_xrays(data, BBBFdict, bbb_obj, BB, models):

    def alpha_OX(log_L2kev):
        """ Relation between accretion disk intrinsic luminosity at 2500 Angstrom and X-rays at 2 keV ."""
        """Lusso&Risaliti +16 gives beta=[0.6-0.65], gamma=[7-8]"""
        beta= 0.643  
        gamma= 6.8734 
        log_L2500A_alphaox= (log_L2kev-gamma)/beta

        return log_L2500A_alphaox

    if models.settings['BBB']=='R06':
        all_bbb_nus, bbb_Fnus_dered = BBBFdict['0.0']                                                    #Intrinsic fluxes without reddening
    elif models.settings['BBB']=='SN12':
        all_bbb_nus, bbb_Fnus_dered = BBBFdict[tuple(np.append(bbb_obj.matched_parkeys_grid[:-1], 0.0))] #Intrinsic fluxes without reddening

    bbb_flux_dered= bbb_Fnus_dered* 10**(BB)
    if BB !=0:
        bbb_flux_dered = bbb_flux_dered*(4*pi*(data.dlum)**2)   ##BB normalization to have units [erg s⁻¹Hz⁻¹] and ranges of values

    """Calculate 2kev (10**17.684 Hz) and 2500 Angstrom (10**15.06 Hz) magnitude in the data and model"""

    flux_2kev = data.fluxes[( 17.284 < (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 18.084 )] ##Rest frame

    if len(flux_2kev) == 0:       #If there the 2keV flux isn't available, use the available flux and assume a power law to estimate it
        flux_Xray = data.fluxes[( 17.60 < (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 19.60 )][0] ##Rest frame
        nu_Xray = data.nus[data.fluxes == flux_Xray][0] 
        nu_2kev = 4.83598*1e17  
        flux_2kev = flux_Xray*(nu_2kev/10**nu_Xray)**(-1.8+1)*np.e**((10**nu_Xray-nu_2kev)/(7.2540*1e19))

    bbb_flux_dered_2500Angs = bbb_flux_dered[(15.04 < all_bbb_nus) & (all_bbb_nus < 15.15 )][0]
    lumfactor = (4. * pi * data.dlum**2.)
    log_L2500A_data_dered = np.log10(bbb_flux_dered_2500Angs)       #Flux in 2500A from data
    log_L2kev_data = np.log10(lumfactor*flux_2kev)                  #Flux in 2keV from data

    log_L2500A_model= alpha_OX(log_L2kev_data)                      #Flux in 2500A from model
    ratio_alpha0x_data= log_L2500A_data_dered - log_L2500A_model

    """Define prior"""
    mu= 0
    sigma= 0.4
    prior_Xrays= Gaussian_prior(mu, sigma, ratio_alpha0x_data) #Promotes accretion disk models that decrease the difference between data and model at 2500A

    return prior_Xrays


def prior_IR_XRays(data, TORUSFdict, tor_obj, TO, models):

    tor_nus, tor_Fnu = TORUSFdict[tor_obj.matched_parkeys_grid]
    tor_flux= tor_Fnu* 10**(TO)
    tor_flux_6microns = tor_flux[(13.59897 < tor_nus) & (tor_nus < 13.79897)][0]  #Flux at 6 microns = 13.69897 log(Hz)
    lumfactor = (4. * pi * data.dlum**2.)
    nuLnu_6microns = (10**13.69897)* tor_flux_6microns * lumfactor                #nuLnu at 6 microns
    x = np.log10(nuLnu_6microns/1e41)
    log_L2_10keV_model = 40.981 + 1.024*x - 0.047*x**2                            #midIR-Xray correlation by Stern 2015 (erg/s)
    logf_2_10keV_model = 22.9494264 + 1.024*x - 0.047*x**2                        #monocromatic flux at 10**17.906 Hz (erg/s/Hz)

    #Central frequency at the band 2-10keV (10**17.906) in rest-frame
    flux2_10keV_data = data.fluxes[( 17.685 < (data.nus+np.log10(1+data.z))) & ((data.nus+np.log10(1+data.z)) < 18.384 )][0]  
    logf2_10keV_data = np.log10(flux2_10keV_data * lumfactor)                    #monocromatic flux at band 2-10keV from data

    ratio_midIR_Xrays= logf2_10keV_data - logf_2_10keV_model

    """Define prior"""
    mu= 0
    sigma= 0.5  
    prior_midIR_Xrays= Gaussian_prior(mu, sigma, ratio_midIR_Xrays)  #Promotes torus models that decrease the difference between data and model at 2-10keV

    return prior_midIR_Xrays


def prior_midIR_UV(data, BBBFdict, bbb_obj, BB, TORUSFdict, tor_obj, TO, models): 

    tor_nus, tor_Fnu = TORUSFdict[tor_obj.matched_parkeys_grid]
    tor_flux= tor_Fnu* 10**(TO)
    tor_flux_6microns = tor_flux[(13.59897 < tor_nus) & (tor_nus < 13.79897)][0]  #Flux at 6 microns = 13.69897 log(Hz)
    lumfactor = (4. * pi * data.dlum**2.)
    x = np.log10(tor_flux_6microns * lumfactor) -27.30103
    log_L2500A_tomodel = (16.2530786 + 1.024*x - 0.047*x**2)/0.643                 #correlations by Stern 2015 + Just et al. 2007

    if models.settings['BBB']=='R06' and models.settings['XRAYS'] != True:
        all_bbb_nus, bbb_Fnus_dered = BBBFdict['0.0']                                                    #Intrinsic fluxes without reddening

    else: 
        EBVbbb_pos = bbb_obj.par_names.index('EBVbbb')
        params = bbb_obj.matched_parkeys_grid
        params[EBVbbb_pos] = str(0.0)
        all_bbb_nus, bbb_Fnus_dered = BBBFdict[tuple(params)]                                            #Intrinsic fluxes without reddening

    bbb_flux_dered= bbb_Fnus_dered* 10**(BB)
    if BB !=0:
        bbb_flux_dered = bbb_flux_dered*(4*pi*(data.dlum)**2)   ##BB normalization to have units [erg s⁻¹Hz⁻¹] and ranges of values
    log_L2500A_bbmodel = np.log10(bbb_flux_dered[( 14.8 < all_bbb_nus) & (all_bbb_nus < 15.15 )][0])        #Flux in 2500A from BB model #15.04

    ratio_2500A= log_L2500A_bbmodel - log_L2500A_tomodel

    """Define prior"""
    mu= 0
    sigma= 0.6 #0.5 (scatter midIR-Xray) + 0.1 (alpha OX)  
    prior_midIR_UV= Gaussian_prior(mu, sigma, ratio_2500A)  #Promotes torus and accretion disk models consistent to each other

    return prior_midIR_UV


def prior_low_AGNfraction(data, models, P, *pars):

    MD = models.dict_modelfluxes
    gal_obj,_,_, bbb_obj = models.dictkey_arrays

    if models.settings['BBB']=='R06': 
        if models.settings['RADIO'] == True:
            GA, SB, TO, BB, RAD = pars[-5:]
        else:
            GA, SB, TO, BB= pars[-4:]
    else:
        if models.settings['RADIO'] == True:
            GA, SB, TO, RAD = pars[-4:]
        else:
            GA, SB, TO = pars[-3:]

    bands, gal_Fnu= MD.GALAXYFdict[gal_obj.matched_parkeys]
    bands, bbb_Fnu = MD.BBBFdict[bbb_obj.matched_parkeys] 

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
    H_0 = 70  #km/s/Mpc
    H_sec = H_0 / 3.0857e19 
    secondsinyear = 31556926

    # Equation for the time elapsed since z and now

    a = 1/(1+z)
    E = O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L
    integrand = lambda z : a/ sqrt(E)        

    #Integration
    z_obs = z
    z_cmb = 1089  

    integral, error = scipy.integrate.quad( integrand , z_obs, z_cmb) #
    
    t = (integral * (1 / H_sec) / secondsinyear) - 1e9  #No stars were born during the first Gyr

    return t
