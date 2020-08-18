

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
import pickle
from astropy.table import Table
from astropy.io import fits, ascii
import scipy
import astropy.constants as const
import astropy.units as u
import itertools
import functions.DICTIONARIES_AGNfitter as dicts
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM       #Cosmology that we assume for estimate luminosity distance                     

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

- Other changes to be done manually (unfortunately):
  (Due to computational time reasons not possible yet to do this automatically)


    (1) Go to fct ymodel in PLOTandWRITE_AGNfitter.py and set (eg. for the galaxy model)
        gal_obj.pick_1D (if your model has one parameter)
        or 
        gal_obj.pick_nD (if your model has more parameters).

    (2) In the same function if your model has more than one parameter, change (e.g. for GALAXYFdict):
        GALAXYFdict[gal_obj.matched_parkeys] (model of one parameter)
        or
        GALAXYFdict[tuple(gal_obj.matched_parkeys)] (model of more than one parameters
"""
def GALAXYfunctions():
    return 0                                                    

def GALAXY(path, modelsettings):
    model_functions = []
    if modelsettings['GALAXY']=='BC03':

        GALAXYFdict_4plot = dict()
        GALAXY_SFRdict = dict()
        GALAXYatt_dict = dict()

        ## Call object containing all galaxy models     
        BC03dict = pickle.load(open(path + 'models/GALAXY/BC03_840seds.pickle', 'rb'), encoding='latin1')    

        ## specify the sizes of the array of parameter values: Here two parameters
        tau_array = BC03dict['tau-values']
        age_array = BC03dict['age-values']
        ebvgal_array = np.array(np.arange(0.,100.,2.5)/100)

        ## produce all combinations of parameter values (indices)
        _, ageidx, tauidx, _, _,_ =  np.shape(BC03dict['SED'])
        
        idxs = [np.arange(ageidx), np.arange(tauidx), np.arange(len(ebvgal_array))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))


        random_array=np.array([int(r) for r in (np.random.rand(42)*16800.)])
        for c in par_idxs_combinations:
                agei=c[0]
                taui=c[1]
                ebvi=c[2]
                #print agei, taui, ebvi
                gal_wl, gal_Fwl =  BC03dict['wavelength'],BC03dict['SED'][:,agei,taui,:,:,:].squeeze()
                gal_nus= gal_wl.to(u.Hz, equivalencies=u.spectral())[::-1] #invert
                gal_Fnu= (gal_Fwl * 3.34e-19 * gal_wl**2.)[::-1]  # Fnu2Flambda
                gal_SFR= BC03dict['SFR'][:,agei,taui,:,:].squeeze()
                gal_nu, gal_Fnu_red = GALAXYred_Calzetti(gal_nus.value[0:len(gal_nus):3], gal_Fnu.value[0:len(gal_nus):3], ebvgal_array[ebvi])  #erg/s/Hz                  
                GALAXYFdict_4plot[str(tau_array.value[taui]),str(np.log10(age_array.value[agei])), str(ebvgal_array[ebvi])] = \
                                                                                        np.log10(gal_nu), renorm_template('GA',gal_Fnu_red)  
                GALAXY_SFRdict[str(tau_array.value[taui]),str(np.log10(age_array.value[agei]))] = gal_SFR 
                gal_Fnu_int = scipy.integrate.trapz(gal_Fnu.value[0:len(gal_nus):3], x=gal_nus.value[0:len(gal_nus):3])
                gal_Fnured_int = scipy.integrate.trapz(gal_Fnu_red, x=gal_nu)
                gal_att_int = gal_Fnu_int- gal_Fnured_int
                GALAXYatt_dict[str(tau_array.value[taui]),str(np.log10(age_array.value[agei])), str(ebvgal_array[ebvi])] = gal_att_int 

        # ## Name the parameters that compose the keyes of the dictionary: GALAXYFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['tau', 'age','EBVgal']
        parameters_types =['grid', 'grid','grid']

        return GALAXYFdict_4plot, GALAXY_SFRdict, GALAXYatt_dict, parameters_names, parameters_types, model_functions


    elif modelsettings['GALAXY']=='BC03_metal':


        GALAXYFdict_4plot = dict()
        GALAXY_SFRdict = dict()
        GALAXYatt_dict = dict()
        ## Call object containing all galaxy models     

        BC03dict = pickle.load(open(path + 'models/GALAXY/BC03_seds_metal_medium.pickle', 'rb'), encoding='latin1')    

        ## specify the sizes of the array of parameter values: Here two parameters
        tau_array = BC03dict['tau-values']
        age_array = BC03dict['age-values']
        metal_array = BC03dict['metallicity-values']
        #ebvgal_array = np.array(np.arange(0.,200.,10.)/100)
        ebvgal_array = np.array(np.arange(0.,60.,5.)/100)

        ## produce all combinations of parameter values (indices)
        metalidx, ageidx, tauidx, _, _,_ =  np.shape(BC03dict['SED'])
        idxs = [np.arange(metalidx), np.arange(ageidx), np.arange(tauidx), np.arange(len(ebvgal_array))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                metali=c[0]
                agei=c[1]
                taui=c[2]
                ebvi=c[3]
                gal_wl, gal_Fwl =  BC03dict['wavelength'],BC03dict['SED'][metali,agei,taui,:,:,:].squeeze()
                gal_nus= gal_wl.to(u.Hz, equivalencies=u.spectral())[::-1]#invert
                gal_Fnu= (gal_Fwl.value * 3.34e-19 * gal_wl**2.)[::-1]  
                gal_nu, gal_Fnu_red = GALAXYred_Calzetti(gal_nus.value[0:len(gal_nus):3], gal_Fnu.value[0:len(gal_nus):3], ebvgal_array[ebvi])                    
                ###!!! gal_Fnu_red
                GALAXYFdict_4plot[str(metal_array[metali]),str(tau_array.value[taui]),str(np.log10(age_array.value[agei])), str(ebvgal_array[ebvi])] = \
                                                                                        np.log10(gal_nu), renorm_template('GA',gal_Fnu_red)       

                gal_SFR= BC03dict['SFR'][metali,agei,taui,:,:].squeeze()
                GALAXY_SFRdict[str(metal_array[metali]),str(tau_array.value[taui]),str(np.log10(age_array.value[agei]))] = gal_SFR         
                gal_Fnu_int = scipy.integrate.trapz(gal_Fnu[0:len(gal_nus):3]*3.826e33, x=gal_nu)
                gal_Fnured_int = scipy.integrate.trapz(gal_Fnu_red*3.826e33, x=gal_nu)
                gal_att_int = gal_Fnu_int.value - gal_Fnured_int
                GALAXYatt_dict[str(metal_array[metali]),str(tau_array.value[taui]),str(np.log10(age_array.value[agei])), str(ebvgal_array[ebvi])] = gal_att_int

        ## Name the parameters that compose the keys of the dictionary: GALAXYFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['metal','tau', 'age','EBVgal']
        parameters_types =['grid','grid', 'grid','grid']

        return  GALAXYFdict_4plot, GALAXY_SFRdict, GALAXYatt_dict, parameters_names, parameters_types, model_functions

def STARBURSTfunctions():
    return 0

def STARBURST(path, modelsettings):
    model_functions = []

    if modelsettings['STARBURST']=='DH02_CE01':

        STARBURSTFdict_4plot = dict()
        STARBURST_LIRdict = dict()

        #Call object containing all starburst models     
        DH02CE01dict = pickle.load(open(path + 'models/STARBURST/DH02_CE01.pickle', 'rb'), encoding='latin1') 
        irlumidx = len(DH02CE01dict['SED'])

        #Construct dictionaries 
        for irlumi in range(irlumidx):
            sb_nu0, sb_Fnu0 = DH02CE01dict['wavelength'][irlumi], DH02CE01dict['SED'][irlumi].squeeze()
            STARBURSTFdict_4plot[str(DH02CE01dict['irlum-values'][irlumi])] = sb_nu0, renorm_template('SB',sb_Fnu0)
            STARBURST_LIRdict[str(DH02CE01dict['irlum-values'][irlumi])] = pow(10,DH02CE01dict['irlum-values'][irlumi])*3.826e33*1e-6

        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['irlum']
        parameters_types =['grid']

        return STARBURSTFdict_4plot, STARBURST_LIRdict, parameters_names, parameters_types, model_functions

    elif modelsettings['STARBURST']=='S17':

        STARBURSTFdict_4plot = dict()
        STARBURST_LIRdict = dict()

        #Call object containing all starburst models     
        dusttable = Table.read(path + 'models/STARBURST/s17_dust.fits')
        pahstable = Table.read(path + 'models/STARBURST/s17_pah.fits')
        
        Dwl, DnuLnu = dusttable['LAM'],dusttable['SED'] #micron, Lsun
        Pwl, PnuLnu = pahstable['LAM'],pahstable['SED'] #micron, Lsun
        Tdust = np.array(dusttable['TDUST'])[0] #K
        fracPAH = np.arange(0.25, 6.25, 0.25)/100
        idxs=[np.arange(len(Tdust)), np.arange(len(fracPAH))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        Dnu= (Dwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        Pnu= (Pwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        DLnu= np.array(DnuLnu[0])/Dnu ###!!!*1e-6 #* u.Lsun.to(u.W)
        PLnu=np.array(PnuLnu[0])/Pnu ###!!!*1e-6 #* u.Lsun.to(u.W)


        #Construct dictionaries 
        for c in par_idxs_combinations:
            t=c[0]
            fp=c[1]

            sb_nu0 = np.array(Dnu[t,:])[::-1]
            sb_Fnu0 = np.array( (1-fracPAH[fp]) * DLnu[t,:] + (fracPAH[fp]) * PLnu[t,:])[::-1]

            STARBURSTFdict_4plot[str(Tdust[t]), str(fracPAH[fp])] = np.log10(sb_nu0), renorm_template('SB',sb_Fnu0)
            STARBURST_LIRdict[str(Tdust[t]), str(fracPAH[fp])] = Lir[t]

        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['Tdust', 'fracPAH']
        parameters_types =['grid', 'grid']

        return STARBURSTFdict_4plot, STARBURST_LIRdict, parameters_names, parameters_types, model_functions

    elif modelsettings['STARBURST']=='S17_newmodel':

        STARBURSTFdict_4plot = dict()
        STARBURST_LIRdict = dict()
        #Call object containing all starburst models     
        dusttable = Table.read(path + 'models/STARBURST/s17_lowvsg_dust.fits')
        pahstable = Table.read(path + 'models/STARBURST/s17_lowvsg_pah.fits')
        
        Dwl, DnuLnu = dusttable['LAM'],dusttable['SED'] #micron, Lsun
        Pwl, PnuLnu = pahstable['LAM'],pahstable['SED'] #micron, Lsun
        Tdust = np.array(dusttable['TDUST'])[0] #K
        Lir=  np.array(dusttable['LIR'])[0] *3.826e33 ###!!!*1e-6 #Lsun2ergs ### consider taking away renormalizaion 1e-6
        fracPAH = np.concatenate(((np.arange(0.0, 0.1, 0.01)/100.),(np.arange(0.1, 5.5, 0.1)/100.)))

        idxs=[np.arange(len(Tdust)), np.arange(len(fracPAH))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        Dnu= (Dwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        Pnu= (Pwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        DLnu= np.array(DnuLnu[0])/Dnu ###!!!*1e-6 #* u.Lsun.to(u.W)
        PLnu=np.array(PnuLnu[0])/Pnu ###!!!*1e-6#* u.Lsun.to(u.W)


        #Construct dictionaries 
        for c in par_idxs_combinations:
            t=c[0]
            fp=c[1]
            #print fracPAH[0]

            sb_nu0 = np.array(Dnu[t,:])[::-1]
            sb_Fnu0 = np.array( (1-fracPAH[fp]) * DLnu[t,:] + (fracPAH[fp]) * PLnu[t,:])[::-1]

            STARBURSTFdict_4plot[str(Tdust[t]), str(fracPAH[fp])] = np.log10(sb_nu0), renorm_template('SB',sb_Fnu0)
            STARBURST_LIRdict[str(Tdust[t]), str(fracPAH[fp])] = Lir[t]
        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['Tdust', 'fracPAH']
        parameters_types =['grid', 'grid']

        return STARBURSTFdict_4plot, STARBURST_LIRdict, parameters_names, parameters_types, model_functions

    elif modelsettings['STARBURST']=='S17_newmodel_radio':

        STARBURSTFdict_4plot = dict()
        STARBURST_LIRdict = dict()

        #Call object containing all starburst models     
        dusttable = Table.read(path + 'models/STARBURST/s17_lowvsg_dust.fits')
        pahstable = Table.read(path + 'models/STARBURST/s17_lowvsg_pah.fits')
         
        Dwl, DnuLnu = dusttable['LAM'],dusttable['SED'] #micron, Lsun
        Pwl, PnuLnu = pahstable['LAM'],pahstable['SED'] #micron, Lsun
        Tdust = np.array(dusttable['TDUST'])[0] #K
        LIR=  np.array(dusttable['LIR'])[0]*3.826e33

        fracPAH = np.concatenate(((np.arange(0.0, 0.1, 0.02)/100.),(np.arange(0.1, 5.5, 0.2)/100.)))
        RADexc= np.logspace(-1, 3, 10) #np.arange(-1.5, 0.9, 0.4)  #np.arange(0, 100, 5) / np.logspace(-1, 3, 10) /np.arange(-0.5, 0.26, 0.152) 

        idxs=[np.arange(len(Tdust)), np.arange(len(fracPAH)),np.arange(len(RADexc))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))
        conv_factor= 3.826e-33 #1e-6
        Dnu= (Dwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        Pnu= (Pwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        DLnu= np.array(DnuLnu[0])/Dnu #*conv_factor #* u.Lsun.to(u.W)
        PLnu=np.array(PnuLnu[0])/Pnu #*conv_factor#* u.Lsun.to(u.W)


        #Construct dictionaries 
        for c in par_idxs_combinations:
            t=c[0]
            fp=c[1]
            re=c[2]
            #print fracPAH[0]

            sb_nu0 = np.array(Dnu[t,:])[::-1]
            sb_Fnu0 = np.array( (1-fracPAH[fp]) * DLnu[t,:] + (fracPAH[fp]) * PLnu[t,:])[::-1]

            rad_sb_nu0 ,rad_sb_Fnu0= RADIO(modelsettings, LIR[t], conv_factor, sb_nu0, sb_Fnu0, RADexc[re]) 

            STARBURSTFdict_4plot[str(Tdust[t]), str(fracPAH[fp]), str(RADexc[re])] = rad_sb_nu0, renorm_template('SB',rad_sb_Fnu0) 
            STARBURST_LIRdict[str(Tdust[t]), str(fracPAH[fp])] = LIR[t]
        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['Tdust', 'fracPAH', 'RADexc']
        parameters_types =['grid', 'grid', 'free']

        return STARBURSTFdict_4plot, STARBURST_LIRdict, parameters_names, parameters_types, model_functions

def BBBfunctions():
    def apply_reddening (bbb_nu, bbb_Fnu, EBV_bbb):
        bbb_nu0, bbb_Fnu_red = BBBred_Prevot(bbb_nu, bbb_Fnu.flatten(), float(EBV_bbb))  
        return bbb_nu0, bbb_Fnu_red
    def add_xrays (bbb_nu, bbb_Fnu, alpha_scat):
        xray_nu, xray_Fnu = XRAYS(modelsettings, bbb_nu, bbb_Fnu.flatten(), float(alpha_scat)) 
        # R06 SED is extended and we need avoid overlapping between BBB template and X-Rays power-law
        bbb_nu0x, bbb_Fnu_redx = np.concatenate((bbb_nu0[bbb_nu0 < 17], xray_nu)), np.concatenate((bbb_Fnu_red[bbb_nu0 < 17], xray_Fnu))
        return bbb_nu0x, bbb_Fnu_redx
    def f0 (): #Dummy function
        return None
    return [apply_reddening,add_xrays,f0]

def BBB(path, modelsettings):
    ## Model from Richards 2006    

    if modelsettings['BBB']=='R06':


        BBBFdict_4plot = dict()
        R06dict = pickle.load(open(path + 'models/BBB/R06.pickle', 'rb'), encoding='latin1') 
        bbb_nu, bbb_Fnu = R06dict['wavelength'], R06dict['SED'].squeeze()
        BBB_functions = BBBfunctions()

        if modelsettings['XRAYS']==False:

            parameters_names =['EBVbbb']
            parameters_types =['free']

            if parameters_types ==['grid']:
                ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)
                model_functions = [0]
                for EBV_bbb in ebvbbb_array:
                    bbb_nu0, bbb_Fnu_red = BBB_functions[0](bbb_nu, bbb_Fnu, EBV_bbb)               
                    BBBFdict_4plot[str(EBV_bbb)] = bbb_nu0, renorm_template('BB', bbb_Fnu_red)
            elif parameters_types ==['free']:
                model_functions = [0]
                ebvbbb_array = np.array([0.,1.0]) # write limits for the free parameter  
                for EBV_bbb in ebvbbb_array:   
                    bbb_nu0, bbb_Fnu_red = BBB_functions[0](bbb_nu, bbb_Fnu, ebvbbb_array[0])  
                    BBBFdict_4plot[str(EBV_bbb)] = bbb_nu0, renorm_template('BB', bbb_Fnu_red)
                
            return BBBFdict_4plot, parameters_names, parameters_types, model_functions

        elif modelsettings['XRAYS']==True:

            parameters_names =['EBVbbb', 'alphaScat'] 
            parameters_types ==['grid', 'grid']

            if parameters_types[0] ==['grid']:
                if parameters_types[1] ==['grid']:
                    model_functions = [2, 1]
                    alpha_scat =  np.arange(-0.2, 0.55, 0.1)
                elif parameters_types[1] ==['free',]:
                    alpha_scat =  np.array([-0.2, 0.55])
                idxs = [np.arange(len(ebvbbb_array)), np.arange(len(alpha_scat))]
                par_idxs_combinations = np.array(list(itertools.product(*idxs)))
                for c in par_idxs_combinations:
                    ebvi=c[0]
                    alpha_scati=c[1]
                    bbb_nu0, bbb_Fnu_red = BBB_functions[0](bbb_nu, bbb_Fnu, ebvi)      
                    if parameters_types[1] ==['grid']:
                        bbb_nu0x, bbb_Fnu_redx = BBB_functions[1](bbb_nu, bbb_Fnu, alpha_scati)
                        BBBFdict_4plot[str(EBV_bbb[ebvi]), str(alpha_scat[alpha_scati])] = bbb_nu0x, renorm_template('BB', bbb_Fnu_redx)      
                    elif parameters_types[1] ==['free']:
                        bbb_nu0x, bbb_Fnu_redx = (bbb_nu0, bbb_Fnu_red)
                        BBBFdict_4plot[str(EBV_bbb[ebvi]), str(alpha_scat[alpha_scati])] = bbb_nu0x, renorm_template('BB', bbb_Fnu_redx)      
            elif parameters_types[0] ==['free']:
                if parameters_types[1] ==['free']:
                    model_functions = [0,1]
                    alpha_scat =  np.array([-0.2, 0.55])
                    ebvbbb_array = np.array([0.,1.0])
                    idxs = [np.arange(len(ebvbbb_array)), np.arange(len(alpha_scat))]
                    par_idxs_combinations = np.array(list(itertools.product(*idxs)))
                    for c in par_idxs_combinations:
                        ebvi=c[0]
                        alpha_scati=c[1]
                        BBBFdict_4plot[str(EBV_bbb[ebvi]), str(alpha_scat[alpha_scati])] = bbb_nu, renorm_template('BB', bbb_Fnu)      

            return BBBFdict_4plot, parameters_names, parameters_types, model_functions

        #Construct dictionaries



    ## Name the parameters that compose the keys of the dictionary: BBFdict_4plot[key]. 
    ## Add the names in the same order as their values are arranged in the dictionary key above.    

    elif modelsettings['BBB']=='SN12':

        BBBFdict_4plot = dict()
        ## Call file containing all galaxy models     
        SN12dict = pickle.load(open(path + 'models/BBB/SN12.pickle', 'rb'), encoding='latin1')    
        parameters_names =['logBHmass', 'logEddra', 'EBVbbb']
        parameters_types =['grid', 'grid', 'grid']


        ## specify the sizes of the array of parameter values: Here two parameters
        ## spin = 0. --> If wished otherwise, request a new modelfile in Github.
        Mbh_array = SN12dict['logBHmass-values']
        EddR_array = SN12dict['logEddra-values']
        ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)
        _, Mbhidx, EddRidx =  np.shape(SN12dict['SED'])

        if modelsettings['XRAYS']==True:
            alpha_scat =  np.arange(-0.2, 0.55, 0.1) 
            parameters_names =['logBHmass', 'logEddra', 'EBVbbb', 'alphaScat']   
    
            ## produce all combinations of parameter values (indices)
            idxs = [np.arange(Mbhidx), np.arange(EddRidx), np.arange(len(ebvbbb_array)), np.arange(len(alpha_scat))]
            par_idxs_combinations = np.array(list(itertools.product(*idxs)))   

            for c in par_idxs_combinations:
                    Mbhi=c[0]
                    EddRi=c[1]
                    ebvi=c[2]
                    alpha_scati = c[3]

                    bbb_nu, bbb_Fnu_nored =  np.log10(SN12dict['frequency']),SN12dict['SED'][:,Mbhi,EddRi].squeeze()
                    bbb_nu0, bbb_Fnu_red = BBBred_Prevot(bbb_nu, bbb_Fnu_nored, ebvbbb_array[ebvi])  

                    #20 differents values for the slope between UV-2keV (scatter)
                    xray_nu, xray_Fnu = XRAYS(modelsettings, bbb_nu, bbb_Fnu_nored, alpha_scat[alpha_scati])
                    bbb_nu0x, bbb_Fnu_redx = np.concatenate((bbb_nu0, xray_nu)), np.concatenate((bbb_Fnu_red, xray_Fnu))

                    BBBFdict_4plot[str(Mbh_array[Mbhi]),str(EddR_array[EddRi]), str(ebvbbb_array[ebvi]), str(alpha_scat[alpha_scati])] = bbb_nu0x, bbb_Fnu_redx

        else:
            parameters_names =['logBHmass', 'logEddra', 'EBVbbb']

            ## produce all combinations of parameter values (indices)
            idxs = [np.arange(Mbhidx), np.arange(EddRidx), np.arange(len(ebvbbb_array))]
            par_idxs_combinations = np.array(list(itertools.product(*idxs)))

            for c in par_idxs_combinations:
                    Mbhi=c[0]
                    EddRi=c[1]
                    ebvi=c[2]

                    bbb_nu, bbb_Fnu_nored =  np.log10(SN12dict['frequency']),SN12dict['SED'][:,Mbhi,EddRi].squeeze()
                    bbb_nu0, bbb_Fnu_red = BBBred_Prevot(bbb_nu, bbb_Fnu_nored, ebvbbb_array[ebvi])  
                    BBBFdict_4plot[str(Mbh_array[Mbhi]),str(EddR_array[EddRi]), str(ebvbbb_array[ebvi])] = bbb_nu0, bbb_Fnu_red
            	      
        return BBBFdict_4plot, parameters_names, parameters_types, model_functions


    elif modelsettings['BBB']=='D12_S':

        BBBFdict_4plot = dict()
        ## Call file containing all galaxy models     
        D12dict = pickle.load(open(path + 'models/BBB/D12_S.pickle', 'rb'), encoding='latin1')    
        parameters_names =['logBHmass', 'logEddra']
        parameters_types =['grid', 'grid']

        ## specify the sizes of the array of parameter values: Here two parameters
        ## spin = 0. --> If wished otherwise, request a new modelfile in Github.
        Mbh_array = D12dict['logBHmass-values']
        EddR_array = D12dict['logEddra-values']

        ## produce all combinations of parameter values (indices)
        _, EddRidx,Mbhidx =  np.shape(D12dict['SED'])
        #idxs = [np.arange(Mbhidx), np.arange(EddRidx), np.arange(len(ebvbbb_array))]
        idxs = [np.arange(EddRidx),np.arange(Mbhidx)]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                EddRi=c[0]
                Mbhi=c[1]
                bbb_nu, bbb_Fnu_nored =  D12dict['frequency'],D12dict['SED'][:,EddRi,Mbhi].squeeze()
                
                ### Apply reddening
                #ebvi=c[2]                
                #bbb_nu, bbb_Fnu_red = BBBred_Prevot(bbb_nu, bbb_Fnu_nored, ebvbbb_array[ebvi])                  
                #BBBFdict_4plot[str(Mbh_array[Mbhi]),str(EddR_array[EddRi]), str(ebvbbb_array[ebvi])] = np.log10(bbb_nu), bbb_Fnu_red        
                
                BBBFdict_4plot[str(Mbh_array[Mbhi]),str(EddR_array[EddRi])] = np.log10(bbb_nu), bbb_Fnu_nored        

        return BBBFdict_4plot, parameters_names, parameters_types, model_functions

    elif modelsettings['BBB']=='D12_K':

        BBBFdict_4plot = dict()
        ## Call file containing all galaxy models     
        D12dict = pickle.load(open(path + 'models/BBB/D12_K.pickle', 'rb'), encoding='latin1')    
        parameters_names =['logBHmass', 'logEddra']
        parameters_types =['grid', 'grid']

        ## specify the sizes of the array of parameter values: Here two parameters
        ## spin = 0. --> If wished otherwise, request a new modelfile in Github.
        Mbh_array = D12dict['logBHmass-values']
        EddR_array = D12dict['logEddra-values']
        #ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)

        ## produce all combinations of parameter values (indices)
        _, EddRidx,Mbhidx =  np.shape(D12dict['SED'])
        #idxs = [np.arange(Mbhidx), np.arange(EddRidx), np.arange(len(ebvbbb_array))]
        idxs = [np.arange(EddRidx),np.arange(Mbhidx)]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                EddRi=c[0]
                Mbhi=c[1]
                bbb_nu, bbb_Fnu_nored =  D12dict['frequency'],D12dict['SED'][:,EddRi,Mbhi].squeeze()
                
                ### Apply reddening
                #ebvi=c[2]                
                #bbb_nu, bbb_Fnu_red = BBBred_Prevot(bbb_nu, bbb_Fnu_nored, ebvbbb_array[ebvi])                  
                #BBBFdict_4plot[str(Mbh_array[Mbhi]),str(EddR_array[EddRi]), str(ebvbbb_array[ebvi])] = np.log10(bbb_nu), bbb_Fnu_red        
                
                BBBFdict_4plot[str(Mbh_array[Mbhi]),str(EddR_array[EddRi])] = np.log10(bbb_nu), bbb_Fnu_nored        

        return BBBFdict_4plot, parameters_names, parameters_types, model_functions

    else:
        print (' ')
        print ('ERROR: The model with the name "'+modelsettings['BBB']+'" does not exist.')

def TORUSfunctions():
    return 0
def TORUS(path, modelsettings):

    model_functions = []
    if modelsettings['TORUS']=='S04':    

        TORUSFdict_4plot  = dict()
        #Call object containing all torus models     
        S04dict = pickle.load(open(path + 'models/TORUS/S04.pickle', 'rb'), encoding='latin1') 
        nhidx=len(S04dict['SED'])
        #Construct dictionaries 
        for nhi in range(nhidx):
            tor_nu0, tor_Fnu0 = S04dict['wavelength'][nhi], S04dict['SED'][nhi].squeeze()
            TORUSFdict_4plot[str(S04dict['Nh-values'][nhi])] = tor_nu0, renorm_template('TO',tor_Fnu0)

        ## Name the parameters that compose the keys of the dictionary: TORUSFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.   
        parameters_names = ['Nh']
        parameters_types = ['grid']

        return TORUSFdict_4plot, parameters_names, parameters_types, model_functions

    elif modelsettings['TORUS']=='NK0': 
        
        TORUSFdict_4plot  = dict()

        NK0dict = pickle.load(open(path + 'models/TORUS/nenkova_v0.pickle', 'rb'), encoding='latin1')  
        incl_idx=len(NK0dict['SED']) 
        #Construct dictionaries 
        for incl_i in range(incl_idx): 

            tor_nu0, tor_Fnu0 = NK0dict['wavelength'][incl_i], NK0dict['SED'][incl_i].squeeze() 
            TORUSFdict_4plot[str(NK0dict['incl-values'][incl_i])] = tor_nu0, renorm_template('TO',tor_Fnu0) 

        parameters_names = ['incl']
        parameters_types = ['grid']

        return TORUSFdict_4plot, parameters_names, parameters_types , model_functions

    elif modelsettings['TORUS']=='SKIRTOR':  #This model has too many parameters --> The code can't find the parameters
        
        TORUSFdict_4plot  = dict()

        SKIRTORdict = pickle.load(open(path + 'models/TORUS/SKIRTOR.pickle', 'rb'), encoding='latin1')  
        
        tv_array = SKIRTORdict['tv-values'].unique()
        p_array = SKIRTORdict['p-values'].unique()
        q_array = SKIRTORdict['q-values'].unique()
        oa_array = SKIRTORdict['oa-values'].unique()
        r_array = SKIRTORdict['r-values'].unique()
        mcl_array = SKIRTORdict['mcl-values'].unique()
        incl_array = SKIRTORdict['incl-values'].unique()

        ## produce all combinations of parameter values (indices)
        idxs = [tv_array, p_array, q_array, oa_array, r_array, mcl_array, incl_array]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                tvi=c[0]
                pi=c[1]
                qi=c[2]
                oai=c[3]
                ri=c[4]
                mcli=c[5]
                incli=c[6]
                model = SKIRTORdict[(SKIRTORdict['tv-values'] == tvi) & (SKIRTORdict['p-values'] == pi)& (SKIRTORdict['q-values'] == qi)& (SKIRTORdict['oa-values'] == oai) & (SKIRTORdict['r-values'] == ri) & (SKIRTORdict['mcl-values'] == mcli) & (SKIRTORdict['incl-values'] == incli)] 
                tor_nu0, tor_Fnu0 =  model['wavelength'].values.item().to_numpy(), model['SED'].values.item().to_numpy()               
                TORUSFdict_4plot[str(tvi),str(pi), str(qi), str(oai), str(ri), str(mcli), str(incli)] = tor_nu0, renorm_template('TO',tor_Fnu0)  

        parameters_names = ['tv', 'p', 'q', 'oa', 'r', 'mcl', 'incl']
        parameters_types = ['grid', 'grid', 'grid', 'grid', 'grid', 'grid', 'grid']
        return TORUSFdict_4plot, parameters_names, parameters_types, model_functions
 
    elif modelsettings['TORUS']=='SKIRTORC': 
        #SKIRTOR model with the parameter values used in X-CIGALE (Yang, Guang, et al. 2020) and inclination as free parameter
        TORUSFdict_4plot  = dict()

        SKIRTORCdict = pickle.load(open(path + 'models/TORUS/SKIRTOR_CIGALE.pickle', 'rb'), encoding='latin1')  
        incl_array = SKIRTORCdict['incl-values']
        #Construct dictionaries 
        for incl_i in incl_array: 

            tor_nu0, tor_Fnu0 = SKIRTORCdict[SKIRTORCdict['incl-values'] == incl_i]['wavelength'].values.item().to_numpy(), SKIRTORCdict[SKIRTORCdict['incl-values'] == incl_i]['SED'].values.item().to_numpy()
            TORUSFdict_4plot[str(incl_i)] = tor_nu0, renorm_template('TO',tor_Fnu0)

        parameters_names = ['incl']
        return TORUSFdict_4plot, parameters_names, parameters_types, model_functions


    elif modelsettings['TORUS']=='SKIRTORM': 
        # SKIRTOR model with averaged SEDs for each inclination
        TORUSFdict_4plot  = dict()

        SKIRTORMdict = pickle.load(open(path + 'models/TORUS/SKIRTOR_mean.pickle', 'rb'), encoding='latin1')  
        incl_array = SKIRTORMdict['incl-values']
        #Construct dictionaries 
        for incl_i in incl_array: 

            tor_nu0, tor_Fnu0 = SKIRTORMdict[SKIRTORMdict['incl-values'] == incl_i]['wavelength'].values.item().to_numpy(), SKIRTORMdict[SKIRTORMdict['incl-values'] == incl_i]['SED'].values.item().to_numpy()
            TORUSFdict_4plot[str(incl_i)] = tor_nu0, renorm_template('TO',tor_Fnu0)

        parameters_names = ['incl']
        return TORUSFdict_4plot, parameters_names, parameters_types , model_functions


    elif modelsettings['TORUS']=='SKIRTORM_2P':
        # SKIRTOR model with averaged SEDs for each inclination and openning angle
        TORUSFdict_4plot  = dict()

        SKIRTORdict = pickle.load(open(path + 'models/TORUS/SKIRTOR_mean_2p.pickle', 'rb'), encoding='latin1')  
        
        oa_array = SKIRTORdict['oa-values'].unique()
        incl_array = SKIRTORdict['incl-values'].unique()

        ## produce all combinations of parameter values (indices)
        idxs = [oa_array, incl_array]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                oai=c[0]
                incli=c[1]
                model = SKIRTORdict[(SKIRTORdict['oa-values'] == oai) & (SKIRTORdict['incl-values'] == incli)] 
                tor_nu0, tor_Fnu0 =  model['wavelength'].values.item().to_numpy(), model['SED'].values.item().to_numpy()               
                TORUSFdict_4plot[str(oai), str(incli)] = tor_nu0, renorm_template('TO',tor_Fnu0)  

        parameters_names = ['oa', 'incl']
        parameters_types = ['grid', 'grid']
        return TORUSFdict_4plot, parameters_names, parameters_types, model_functions


    elif modelsettings['TORUS']=='SKIRTORM_3P':
        # SKIRTOR model with averaged SEDs for each inclination, openning angle and optical depth
        TORUSFdict_4plot  = dict()

        SKIRTORdict = pickle.load(open(path + 'models/TORUS/SKIRTOR_mean_3p.pickle', 'rb'), encoding='latin1')  
        
        oa_array = SKIRTORdict['oa-values'].unique()
        incl_array = SKIRTORdict['incl-values'].unique()
        tv_array = SKIRTORdict['tv-values'].unique()

        ## produce all combinations of parameter values (indices)
        idxs = [oa_array, incl_array, tv_array]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                oai=c[0]
                incli=c[1]
                tvi = c[2]
                model = SKIRTORdict[(SKIRTORdict['oa-values'] == oai) & (SKIRTORdict['incl-values'] == incli) & (SKIRTORdict['tv-values'] == tvi)] 
                tor_nu0, tor_Fnu0 =  model['wavelength'].values.item().to_numpy(), model['SED'].values.item().to_numpy()               
                TORUSFdict_4plot[str(oai), str(incli), str(tvi)] = tor_nu0, renorm_template('TO',tor_Fnu0)  

        parameters_names = ['oa', 'incl', 'tv']
        parameters_types = ['grid', 'grid', 'grid']
        return TORUSFdict_4plot, parameters_names, parameters_types, model_functions

    elif modelsettings['TORUS']=='SKIRTORM_4P':
        # SKIRTOR model with averaged SEDs for each inclination, openning angle and optical depth
        TORUSFdict_4plot  = dict()

        SKIRTORdict = pickle.load(open(path + 'models/TORUS/SKIRTOR_mean_4p.pickle', 'rb'), encoding='latin1')  
        
        oa_array = SKIRTORdict['oa-values'].unique()
        incl_array = SKIRTORdict['incl-values'].unique()
        tv_array = SKIRTORdict['tv-values'].unique()
        p_array = SKIRTORdict['p-values'].unique()

        ## produce all combinations of parameter values (indices)
        idxs = [oa_array, incl_array, tv_array, p_array]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                oai=c[0]
                incli=c[1]
                tvi = c[2]
                pi = c[3]
                model = SKIRTORdict[(SKIRTORdict['oa-values'] == oai) & (SKIRTORdict['incl-values'] == incli) & (SKIRTORdict['tv-values'] == tvi) & (SKIRTORdict['p-values'] == pi)] 
                tor_nu0, tor_Fnu0 =  model['wavelength'].values.item().to_numpy(), model['SED'].values.item().to_numpy()               
                TORUSFdict_4plot[str(oai), str(incli), str(tvi), str(pi)] = tor_nu0, renorm_template('TO',tor_Fnu0)  

        parameters_names = ['oa', 'incl', 'tv', 'p']
        parameters_types = ['grid', 'grid', 'grid', 'grid']
        return TORUSFdict_4plot, parameters_names, parameters_types, model_functions




def RADIO(modelsettings, LIR, conv_factor, sb_nu0, sb_Fnu0, RAD_excess):

    if modelsettings['RADIO']==True:  

        #q_IR_r14 = np.random.normal(2.64, 0.26,1)
        q_IR_r14 = 2.64 
        alpha_syn  = -0.75
        alpha_th = -0.1
        nth_frac =0.9

        L14 = 10**(np.log10(LIR*conv_factor)-q_IR_r14)/3.75e12#1.4e9 #to Wats
        nu_spacing= (np.log10(sb_nu0)[2]-np.log10(sb_nu0)[1])
        radio_points = (np.log10(sb_nu0)[0]-9)/nu_spacing
        radio_nu14= np.log10(1.4e9)

        radio_nu = np.arange(np.log10(sb_nu0)[0]- nu_spacing*int(radio_points),np.log10(sb_nu0)[0], nu_spacing)
        #radio_nu2 = np.concatenate((radio_nu, np.log10(sb_nu0)[:np.argmax(sb_Fnu0)]))
        radio_SB_nu = np.concatenate((radio_nu, np.log10(sb_nu0)[:np.argmax(sb_Fnu0)]))
        #all_nu=  np.concatenate((radio_nu, np.log10(sb_nu0)))
        optical_nu = np.arange(np.log10(sb_nu0)[-1], 16, nu_spacing)
        radio_opt_nu = np.concatenate((radio_nu, np.log10(sb_nu0), optical_nu))


        #Lsb = np.concatenate((sb_Fnu0[0]*1e-4*np.ones(len(radio_nu)),sb_Fnu0))
        Lsb = np.concatenate((sb_Fnu0[0]*1e-4*np.ones(len(radio_nu)),sb_Fnu0, sb_Fnu0[-1]*1e-4*np.ones(len(optical_nu))))

        Lsyn_0 = 10**(-1*alpha_syn* np.log10(1.4e9/10**radio_SB_nu[0])) * L14*(nth_frac)
        #Lsyn_rad = Lsyn_0 * 10**(alpha_syn* np.log10(10**radio_nu_th/10**radio_nu_th[0])) 
        #Lsyn = np.concatenate((Lsyn_rad, 0*np.ones(len(all_nu)-len(Lsyn_rad))))
        Lsyn_SB = Lsyn_0 * 10**(alpha_syn* np.log10(10**radio_SB_nu/10**radio_SB_nu[0])) 
        Lsyn = np.concatenate((Lsyn_SB, 0*np.ones(len(radio_opt_nu)-len(Lsyn_SB))))

        #Lth_0 = 10**(-1*alpha_th* np.log10(1.4e9/10**radio_nu2[0])) * L14*(1.-nth_frac)
        Lth_0 = 10**(-1*alpha_th* np.log10(1.4e9/10**radio_SB_nu[0])) * L14*(1.-nth_frac)
        Lth_SB= Lth_0* 10**(alpha_th* np.log10(10**radio_SB_nu/10**radio_SB_nu[0]))  
        #Lth_rad= Lth_0* 10**(alpha_th* np.log10(10**radio_nu2/10**radio_nu2[0])) 
        #Lth = np.concatenate((Lth_rad, 0*np.ones(len(all_nu)-len(Lth_rad))))
        Lth = np.concatenate((Lth_SB, 0*np.ones(len(radio_opt_nu)-len(Lth_SB))))

        Lsyn_AGN = RAD_excess*L14*10**(alpha_syn* np.log10(10**radio_opt_nu/10**radio_opt_nu[0]))


        Lir_rad= Lsb+Lsyn+Lth+Lsyn_AGN

        return  radio_opt_nu, Lir_rad #radio_nu_syn, Lir_rad   all_nu, Lir_rad

    else:

        print ('No radio data included in the fit.')


def XRAYS(bbb_nu, bbb_Fnu, scatter):
    
    f = interp1d(bbb_nu, bbb_Fnu, kind = 'nearest', bounds_error=False, fill_value=0.) 
    nu_2500 = (3*1e8)/(2500*1e-10)                               # frequency at 2500 Angstroms
    L_2500 = f(np.log10(nu_2500))                                # Luminosity at 2500 Angstroms
    mean_alpha = -0.137*np.log10(L_2500) + 2.638                 # alpha_OX-L_2500 relation
    alpha = mean_alpha + scatter                                 # Scatter in alpha_OX-L_2500 (-2sigma, 2sigma)

    nu_2kev = 4.83598*1e17                                       # frequency at 2 keV
    Fnu_2kev = L_2500*10**(alpha*(np.log10(nu_2kev/nu_2500)))    # Luminosity at 2keV

    #Proportionality constant a to scale x-ray power-law in 2keV to the value found with alpha_OX-L_2500
    h = 4.135667731*1e-15*1e-3                                   #eV/Hz --> keV/Hz
    Gamma = 1.8
    a = Fnu_2kev/((h*nu_2kev)**(-Gamma+1)*np.e**(-nu_2kev/(7.2540*1e19)))

    xray_nu = np.logspace(17, 19, 1000)                         #with a hole between BB template and X-Rays
    xray_Fnu = a*(h*xray_nu)**(-Gamma+1)*np.e**(-xray_nu/(7.2540*1e19))

    return np.log10(xray_nu), xray_Fnu




"""===================================================
Reddening functions    
==================================================="""

# def REDDENING(modelsettings):

#     if modelsettings['gal_reddening']=='Calzetti': 

#         """
#         This function computes the effect of reddening in the galaxy template (Calzetti law)

#         ## input:
#         -frequencies in log nu
#         - Fluxes in Fnu
#         - the reddening value E(B-V)_gal
#         ## output:

#         """
#         RV = 4.05        

#         c =2.998 * 1e8 
#         gal_lambda_m = c / gal_nu * 1e6#in um 
#         wl = gal_lambda_m[::-1]  #invert for lambda
#         k = np.zeros(len(wl))

#         w0 = tuple([wl <= 0.12])
#         w1 = tuple([wl < 0.63])
#         w2 = tuple([wl >= 0.63])

#         x1 = np.argmin(np.abs(wl - 0.12))
#         x2 = np.argmin(np.abs(wl - 0.125))

#         k[w2] = 2.659 * (-1.857 + 1.040 /wl[w2])+RV
#         k[w1] = 2.659 * (-2.156 + (1.509/wl[w1]) - (0.198/wl[w1]**2) + (0.011/wl[w1]**3))+RV
#         k[w0] = k[x1] + ((wl[w0] - 0.12) * (k[x1] - k[x2]) / (wl[x1] - wl[x2])) +RV


#         gal_k= k[::-1] #invert for nus
#         gal_Fnu_red = gal_Fnu* 10**(-0.4 * gal_k * GAebv)
#         return gal_nu, gal_Fnu_red
#     if modelsettings['bbb_reddening']=='Prevot_SMC': 
#         """
        
#         ## input:

#         ## output:

#         """
#         #Application of reddening - reading E(B-V) from MCMC sampler
#         RV= 2.72

#         #converting freq to wavelength, to be able to use prevots function instead on simple linear interpolation 
#         redd_x =  2.998 * 1e10 / (10**(bbb_x)* 1e-8)
#         redd_x= redd_x[::-1]

#         #    Define prevots function for the reddening law redd_k    
#         def function_prevot(x, RV):
#                y=1.39*pow((pow(10.,-4.)*x),-1.2)-0.38 ;
#                return y 

#         bbb_k = function_prevot(redd_x, RV)

#         bbb_k= bbb_k[::-1]

#         bbb_Lnu_red = bbb_y * 10**(-0.4 * bbb_k * BBebv)

#         bbb_Lnu_red[np.isnan(bbb_Lnu_red)]=bbb_y[np.isnan(bbb_Lnu_red)]

#         return bbb_x, bbb_Lnu_red


def BBBred_Prevot(bbb_x, bbb_y, BBebv ):

    """
    This function computes the effect of reddening in the accretion disk template (Prevot law for Small Magellanic Cloud)

    ## input:
    -frequencies in log nu
    - Fluxes in Fnu
    - the reddening value E(B-V)_bb
    ## output:

    """
    #Application of reddening - reading E(B-V) from MCMC sampler
    RV= 2.72

    #converting freq to wavelength, to be able to use prevots function instead on simple linear interpolation 
    redd_x =  2.998 * 1e10 / (10**(bbb_x)* 1e-8)
    redd_x= redd_x[::-1]

    #    Define prevots function for the reddening law redd_k    
    def function_prevot(x, RV):
           y=1.39*pow((pow(10.,-4.)*x),-1.2)-0.38 ;
           return y 

    bbb_k = function_prevot(redd_x, RV)
    bbb_k= bbb_k[::-1]
    bbb_Lnu_red = bbb_y * 10**(-0.4 * bbb_k * BBebv)
    bbb_Lnu_red[np.isnan(bbb_Lnu_red)]=bbb_y[np.isnan(bbb_Lnu_red)]

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

    w0 = tuple([wl <= 0.12])
    w1 = tuple([wl < 0.63])
    w2 = tuple([wl >= 0.63])

    x1 = np.argmin(np.abs(wl - 0.12))
    x2 = np.argmin(np.abs(wl - 0.125))

    k[w2] = 2.659 * (-1.857 + 1.040 /wl[w2])+RV
    k[w1] = 2.659 * (-2.156 + (1.509/wl[w1]) - (0.198/wl[w1]**2) + (0.011/wl[w1]**3))+RV
    k[w0] = k[x1] + ((wl[w0] - 0.12) * (k[x1] - k[x2]) / (wl[x1] - wl[x2])) +RV


    gal_k= k[::-1] #invert for nus
    gal_Fnu_red = gal_Fnu* 10**(-0.4 * gal_k * GAebv)
    return gal_nu, gal_Fnu_red


def GALAXYred_CharlotFall(gal_nu, gal_Fnu,GAebv):

    """
    This function computes the effect of reddening in the galaxy template (Charlot and Fall +00)
    ## input:
    -frequencies in log nu
    - Fluxes in Fnu
    - the reddening value E(B-V)_gal
    ## output:

    """
    RV = 5.9        

    c =2.998 * 1e8 
    gal_lambda_m = c / gal_nu * 1e6#in um 
    wl = gal_lambda_m[::-1]  #invert for lambda
    k = np.zeros(len(wl))

    kcf = RV * (wl/5500)**(-0.7)

    gal_k= kcf[::-1] #invert for nus
    gal_Fnu_red = gal_Fnu* 10**(-0.4 * gal_k * GAebv)
    return gal_nu, gal_Fnu_red



Angstrom = 1e10

def z2Dlum(z):

    cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.266)   
    dlum_cm = cosmo.luminosity_distance(z).to(u.cm).value
    dlum_Mpc = dlum_cm/3.08567758e24 

    return dlum_cm
   
def fluxlambda_2_fluxnu (flux_lambda, wl_angst):

    c = 2.99792458e8

    flux_nu = flux_lambda * (wl_angst**2. ) / c /Angstrom
    return flux_nu

"""---------------------------------------------
             COMPUTED QUANTITIES
-----------------------------------------------"""

def stellar_info(chain, data, models):

    """
    computes stellar masses and SFRs
    """

    gal_obj,_,_,_ = models.dictkey_arrays
    MD= models.dict_modelfluxes

    if len(gal_obj.par_names)==3:
        tau_mcmc = chain[:,0]  
        age_mcmc = chain[:,1] 
    elif len(gal_obj.par_names)==4:
        metal_mcmc = chain[:,0] 
        tau_mcmc = chain[:,1]     
        age_mcmc = chain[:,2] 
    GA = chain[:, -4] #- 18. ## 1e18 is the common normalization factor used in parspace.ymodel 
                            ## in order to have comparable NORMfactors    
    z = data.z
    distance = z2Dlum(z)

    #constants
    solarlum = const.L_sun.to(u.erg/u.second) #3.839e33
    solarmass = const.M_sun

    Mstar_list=[]
    SFR_list=[]

    for i in range (len (tau_mcmc)):        
        N = 10**GA[i]* 4* pi* distance**2 / (solarlum.value)/ (1+z)
        N = renorm_template('GA', N)

        if len(gal_obj.par_names)==3:
###!        gal_obj.pick_nD(tuple([tau_mcmc[i], age_mcmc[i], 0.]))
            gal_obj.pick_nD(tuple([tau_mcmc[i], age_mcmc[i], 0.]))
            tau_dct, age_dct, ebvg_dct=gal_obj.matched_parkeys
            SFR_mcmc =MD.GALAXY_SFRdict[tau_dct, age_dct]
        elif len(gal_obj.par_names)==4:
###!        gal_obj.pick_nD(tuple([metal_mcmc[i], tau_mcmc[i], age_mcmc[i], 0.]))
            gal_obj.pick_nD(tuple([metal_mcmc[i], tau_mcmc[i], age_mcmc[i], 0.]))
            metal_dct,tau_dct, age_dct, ebvg_dct=gal_obj.matched_parkeys
            SFR_mcmc =MD.GALAXY_SFRdict[metal_dct,tau_dct, age_dct]

        # Calculate Mstar. BC03 templates are normalized to M* = 1 M_sun. 
        # Thanks to Kenneth Duncan, and his python version of BC03, smpy
        Mstar = np.log10(N * 1) 
        #Calculate SFR. output is in [Msun/yr]. 
        #print tau_mcmc[i], age_mcmc[i], (N * SFR_mcmc), np.log10(N), SFR_mcmc 
        SFR = N * SFR_mcmc
        SFR_list.append(SFR.value)    
        Mstar_list.append(Mstar)    

    return np.array(Mstar_list) , np.array(SFR_list)


def stellar_info_array(chain_flat, data, models, Nthin_compute):

    """
    computes arrays of stellar masses and SFRs
    """
    import random
    Ns, Npar = np.shape(chain_flat) 
    chain_thinned = chain_flat[random.sample(list(np.arange(Ns)), Nthin_compute),:]

    Mstar, SFR = stellar_info(chain_thinned, data, models)
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



"""---------------------------------------------
            RENORMALIZATION OF MODELS
-----------------------------------------------"""

def renorm_template(model, Fnu):
    if model=='GA':
        Fnu_norm = Fnu/1e18
        return Fnu_norm
    elif model== 'SB':
        Fnu_norm = Fnu/1e20
        return Fnu_norm
    elif model== 'TO':
        Fnu_norm = Fnu/1e-40
        return Fnu_norm
    elif model== 'BB':
        Fnu_norm = Fnu/1e60 ## 1e60 change to 1e64

        return Fnu_norm


class MODELS: #is this used somewhere else?

    def __init__(self, z, models_settings):
        self.name= 'models'
        self.z=z
        self.settings = models_settings

    def DICTS(self, filters, Modelsdict):
        """
        Helps transporting the dictionary content
        corresponding to the redshift of the source
        """
        # z_array = np.array(list(Modelsdict.keys()))###erase
        # idx = (np.abs(z_array.astype(float)-self.z)).argmin()###erase
        # z_key = z_array[idx]###erase
        #self.dict_modelfluxes = Modelsdict[z_key]
        self.dict_modelfluxes = Modelsdict#[z_key]
        self.dictkey_arrays = dicts.dictkey_arrays(self.dict_modelfluxes)
        self.dictkey_arrays_4plot = dicts.dictkey_arrays_4plot(self.dict_modelfluxes)

