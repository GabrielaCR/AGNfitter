

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
from astropy.io import fits, ascii
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

- Other changes to be done manually (unfortunately):
  (Due to computational time reasons not possible yet to do this automatically)

    (1) Go to fct ymodel in PARAMETERSPACE_AGNfitter.py and set (eg. for the galaxy model)
        gal_obj.pick_1D (if your model has one parameter)
        or 
        gal_obj.pick_nD (if your model has more parameters).

    (2) In the same function if your model has more than one parameter, change (e.g. for GALAXYFdict):
        GALAXYFdict[gal_obj.matched_parkeys] (model of one parameter)
        or
        GALAXYFdict[tuple(gal_obj.matched_parkeys)] (model of more than one parameters

"""
                                                         

def GALAXY(path, modelsettings):

    if modelsettings['GALAXY']=='BC03':


        GALAXYFdict_4plot = dict()
        GALAXY_SFRdict = dict()
        ## Call object containing all galaxy models     
        #BC03dict = cPickle.load(file(path + 'models/GALAXY/BC03_840seds.pickle', 'rb'))    
        BC03dict = cPickle.load(file(path + 'models/GALAXY/BC03_840seds.pickle', 'rb'))    

        ## specify the sizes of the array of parameter values: Here two parameters
        tau_array = BC03dict['tau-values']
        age_array = BC03dict['age-values']
        ebvgal_array = np.array(np.arange(0.,200.,5.)/200)

        ## produce all combinations of parameter values (indices)
        _, ageidx, tauidx, _, _,_ =  np.shape(BC03dict['SED'])
        
        idxs = [np.arange(ageidx), np.arange(tauidx), np.arange(len(ebvgal_array))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        ##########################################erase!!!!!
        i_1= 0 
        i_2=0#erase!!!!!
        SFR_array = BC03dict['SFR']
        import random
        mock_gal_catalog=[]
        ##########################################erase!!!!!

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
                GALAXY_SFRdict[str(tau_array.value[taui]),str(np.log10(age_array.value[agei]))] = gal_SFR         
                gal_nu, gal_Fnu_red = GALAXYred_Calzetti(gal_nus.value[0:len(gal_nus):3], gal_Fnu.value[0:len(gal_nus):3], ebvgal_array[ebvi])  #erg/s/Hz                  
                GALAXYFdict_4plot[str(tau_array.value[taui]),str(np.log10(age_array.value[agei])), str(ebvgal_array[ebvi])] = \
                                                                                        np.log10(gal_nu), renorm_template('GA',gal_Fnu_red)        
                ##########################################erase!!!!!
                i_1=i_1+1  
                if i_1 in random_array:
                    i_2=i_2+1
                    zm=random.randint(1,20)*0.1
                    GAm = random.randint(2,4)
                    dm = z2Dlum(zm)
                    solarlum = const.L_sun.to(u.erg/u.second) #3.839e33
                    Nm = 10**GAm * 4* pi* dm**2 / (solarlum.value)/ (1+zm)
                    Nm =renorm_template('GA', Nm)
                    xm= gal_nu
                    ym= gal_Fnu_red * Nm #*solarlum.value *1e-23 *1e-3
                    SFRm = BC03dict['SFR'][:,agei,taui,:,:][0][0][0].value * Nm
                    #dlum= z2Dlum(1)
                    #lumfactor = (4. * pi * dlum**2.)

                    ###Choose data to simulate

                    ### look for this frequencies
                    #GALEX_1500|1516|15.296
                    #u_SDSS|3591|14.921
                    #g_SDSS|4723|14.802|
                    #r_SDSS|6213|14.683|
                    #i_SDSS|7523|14.600|
                    #z_SDSS|8855|14.529|
                    #F105W|10517|14.455
                    #F140W|13877|14.334
                    #F160W|15348|14.291
                    #H_2mass|16467|14.260
                    #Ks_2mass|21641|14.141
                    #IRAC1|35375|13.928|
                    #IRAC2|44751|13.826|

                    ### Calculate fluxes to extract from model. Redshift model and extract at telescope bands.
                    freq_obs_array=np.array([13.836, 13.928, 14.141, 14.260, 14.291, 14.331, 14.455, 14.529, 14.600, 14.683, 14.802, 14.921, 15.296])
                    freq_intrinsic_array= np.log10(10**freq_obs_array*(1+zm))
                    fidxs=[]
                    for fi in freq_intrinsic_array:        
                        fidx = (np.abs(np.log10(xm)-fi)).argmin()
                        fidxs.append(fidx)
                    fidxsa= np.array(fidxs)[0:2].tolist()
                    fidxsb= np.array(fidxs)[2:].tolist()
                    ymerr=np.hstack((ym[fidxsa]*0.5, ym[fidxsb]*0.1))

                    ### Construct output data, ID, SFR, z, nu, Fnu, Fnu_err
                    data =np.array([np.hstack((i_2,SFRm, zm, 10**freq_obs_array, ym[fidxs], ym[fidxs]*0.1))])

                    ### Do plots to compare
                    fig = plt.figure()
                    plt.errorbar(freq_obs_array, ym[fidxs]*(10**freq_obs_array), yerr= ymerr*(10**freq_obs_array), ecolor='k', linestyle='', ms=3)
                    plt.plot(np.log10(xm/(1+zm)),ym*(xm/(1+zm)),'r-',lw=0.5, label=' age='+str(np.log10(age_array.value[agei])) +\
                                                                                    ' tau='+str(tau_array.value[taui])+ \
                                                                                    ' EB-V='+str(ebvgal_array[ebvi]) +\
                                                                                    ' SFR='+str(SFRm)+', z='+str(zm))
                    plt.yscale('log')
                    plt.legend(prop={'size':7})
                    #plt.savefig('/Users/Gabriela/Desktop/AGNfitter/data/mockGalSED_'+str(i_2)+'.png')
                    plt.close

                    mock_gal_catalog.append(data)
                
        mock_gal_catalog= np.reshape(np.array(mock_gal_catalog),[42,42])

        ### Write mock catalog
        # ascii.write(mock_gal_catalog, '/Users/Gabriela/Desktop/AGNfitter/data/mockGalSED.dat', names=['N','SFR', 'z', \
        #                                                                                                 'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12','x13',\
        #                                                                                                 'y1','y2','y3','y4','y5','y6','y7','y8','y9','y10','y11','y12','y13',\
        #                                                                                                 'ye1','ye2','ye3','ye4','ye5','ye6','ye7','ye8','ye9','ye10','ye11','ye12','ye13'],overwrite=True)
        #########################################

        # ## Name the parameters that compose the keyes of the dictionary: GALAXYFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['tau', 'age','EBVgal']

        #plt.show()


         
        return  GALAXYFdict_4plot, GALAXY_SFRdict, parameters_names


    elif modelsettings['GALAXY']=='BC03_metal':


        GALAXYFdict_4plot = dict()
        GALAXY_SFRdict = dict()
        ## Call object containing all galaxy models     

        BC03dict = cPickle.load(file(path + 'models/GALAXY/BC03_seds_metal_medium.pickle', 'rb'))    

        ## specify the sizes of the array of parameter values: Here two parameters
        tau_array = BC03dict['tau-values']
        age_array = BC03dict['age-values']
        metal_array = BC03dict['metallicity-values']
        ebvgal_array = np.array(np.arange(0.,250.,25.)/100)

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
                gal_Fnu= (gal_Fwl * 3.34e-19 * gal_wl**2.)[::-1]  
                gal_SFR= BC03dict['SFR'][metali,agei,taui,:,:].squeeze()
                GALAXY_SFRdict[str(metal_array[metali]),str(tau_array.value[taui]),str(np.log10(age_array.value[agei]))] = gal_SFR         
                gal_nu, gal_Fnu_red = GALAXYred_Calzetti(gal_nus.value[0:len(gal_nus):3], gal_Fnu.value[0:len(gal_nus):3], ebvgal_array[ebvi])                    
                GALAXYFdict_4plot[str(metal_array[metali]),str(tau_array.value[taui]),str(np.log10(age_array.value[agei])), str(ebvgal_array[ebvi])] = \
                                                                                        np.log10(gal_nu), renorm_template('GA',gal_Fnu_red)       

        ## Name the parameters that compose the keys of the dictionary: GALAXYFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['metal','tau', 'age','EBVgal']

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
            STARBURSTFdict_4plot[str(DH02CE01dict['irlum-values'][irlumi])] = sb_nu0, renorm_template('SB',sb_Fnu0)

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
        fracPAH = np.arange(0.25, 6.25, 0.25)/100
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

            STARBURSTFdict_4plot[str(Tdust[t]), str(fracPAH[fp])] = np.log10(sb_nu0), renorm_template('SB',sb_Fnu0)

        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['Tdust', 'fracPAH']

        return STARBURSTFdict_4plot, parameters_names

    elif modelsettings['STARBURST']=='S17_newmodel':

        STARBURSTFdict_4plot = dict()

        #Call object containing all starburst models     
        dusttable = Table.read(path + 'models/STARBURST/s17_lowvsg_dust.fits')
        pahstable = Table.read(path + 'models/STARBURST/s17_lowvsg_pah.fits')
        
        Dwl, DnuLnu = dusttable['LAM'],dusttable['SED'] #micron, Lsun
        Pwl, PnuLnu = pahstable['LAM'],pahstable['SED'] #micron, Lsun
        Tdust = np.array(dusttable['TDUST'])[0] #K
        #fracPAH = np.arange(0.0, 5.5, 0.05)/100
        fracPAH = np.concatenate(((np.arange(0.0, 0.1, 0.01)/100.),(np.arange(0.1, 5.5, 0.1)/100.)))


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
            #print fracPAH[0]

            sb_nu0 = np.array(Dnu[t,:])[::-1]
            sb_Fnu0 = np.array( (1-fracPAH[fp]) * DLnu[t,:] + (fracPAH[fp]) * PLnu[t,:])[::-1]
            # print type(Tdust)#[t]
            # print (Dnu[t,:])
            # print np.shape(sb_Fnu0), sb_Fnu0

            STARBURSTFdict_4plot[str(Tdust[t]), str(fracPAH[fp])] = np.log10(sb_nu0), renorm_template('SB',sb_Fnu0)

        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['Tdust', 'fracPAH']

        return STARBURSTFdict_4plot, parameters_names

    elif modelsettings['STARBURST']=='S17_newmodel_radio':

        STARBURSTFdict_4plot = dict()

        #Call object containing all starburst models     
        dusttable = Table.read(path + 'models/STARBURST/s17_lowvsg_dust.fits')
        pahstable = Table.read(path + 'models/STARBURST/s17_lowvsg_pah.fits')
         
        Dwl, DnuLnu = dusttable['LAM'],dusttable['SED'] #micron, Lsun
        Pwl, PnuLnu = pahstable['LAM'],pahstable['SED'] #micron, Lsun
        Tdust = np.array(dusttable['TDUST'])[0] #K
        LIR=  np.array(dusttable['LIR'])[0]

        #fracPAH = np.arange(0.0, 5.5, 0.05)/100
        fracPAH = np.concatenate(((np.arange(0.0, 0.1, 0.01)/100.),(np.arange(0.1, 5.5, 0.1)/100.)))
        RADexc= np.arange(0, 100, 5)

        idxs=[np.arange(len(Tdust)), np.arange(len(fracPAH)),np.arange(len(RADexc))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))
        conv_factor=1e-6
        Dnu= (Dwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        Pnu= (Pwl[0] * u.micron).to(u.Hz, equivalencies=u.spectral())
        DLnu= np.array(DnuLnu[0])/Dnu *conv_factor #* u.Lsun.to(u.W)
        PLnu=np.array(PnuLnu[0])/Pnu *conv_factor#* u.Lsun.to(u.W)


        #Construct dictionaries 
        for c in par_idxs_combinations:
            t=c[0]
            fp=c[1]
            re=c[2]
            #print fracPAH[0]

            sb_nu0 = np.array(Dnu[t,:])[::-1]
            sb_Fnu0 = np.array( (1-fracPAH[fp]) * DLnu[t,:] + (fracPAH[fp]) * PLnu[t,:])[::-1]

            rad_sb_nu0 ,rad_sb_Fnu0= RADIO(LIR[t], conv_factor, sb_nu0, sb_Fnu0, RADexc[re])

            STARBURSTFdict_4plot[str(Tdust[t]), str(fracPAH[fp]), str(RADexc[re])] = rad_sb_nu0, renorm_template('SB',rad_sb_Fnu0)

        ## Name the parameters that compose the keys of the dictionary: STARBURSTFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.    
        parameters_names =['Tdust', 'fracPAH', 'RADexc']

        return STARBURSTFdict_4plot, parameters_names

def BBB(path, modelsettings):

    if modelsettings['BBB']=='R06':

        BBBFdict_4plot = dict()
        R06dict = cPickle.load(file(path + 'models/BBB/R06.pickle', 'rb')) 
        parameters_names =['EBVbbb']
        ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)

        bbb_nu, bbb_Fnu = R06dict['wavelength'], R06dict['SED'].squeeze()
        
        #Construct dictionaries
        for EBV_bbb in ebvbbb_array:
            bbb_nu0, bbb_Fnu_red = BBBred_Prevot(bbb_nu, bbb_Fnu, EBV_bbb )
            BBBFdict_4plot[str(EBV_bbb)] = bbb_nu0, renorm_template('BB', bbb_Fnu_red)
        return BBBFdict_4plot, parameters_names


    ## Name the parameters that compose the keys of the dictionary: BBFdict_4plot[key]. 
    ## Add the names in the same order as their values are arranged in the dictionary key above.    

    elif modelsettings['BBB']=='SN12':

        BBBFdict_4plot = dict()
        ## Call file containing all galaxy models     
        SN12dict = cPickle.load(file(path + 'models/BBB/SN12.pickle', 'rb'))    
        parameters_names =['logBHmass', 'logEddra', 'EBVbbb']#SN12dict['parameters']

        ## specify the sizes of the array of parameter values: Here two parameters
        ## spin = 0. --> If wished otherwise, request a new modelfile in Github.
        Mbh_array = SN12dict['logBHmass-values']
        EddR_array = SN12dict['logEddra-values']
        ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)

        ## produce all combinations of parameter values (indices)
        _, Mbhidx, EddRidx =  np.shape(SN12dict['SED'])
        idxs = [np.arange(Mbhidx), np.arange(EddRidx), np.arange(len(ebvbbb_array))]
        par_idxs_combinations = np.array(list(itertools.product(*idxs)))

        for c in par_idxs_combinations:
                Mbhi=c[0]
                EddRi=c[1]
                ebvi=c[2]
                bbb_nu, bbb_Fnu_nored =  SN12dict['frequency'],SN12dict['SED'][:,Mbhi,EddRi].squeeze()
                bbb_nu, bbb_Fnu_red = BBBred_Prevot(bbb_nu, bbb_Fnu_nored, ebvbbb_array[ebvi])                  
                BBBFdict_4plot[str(Mbh_array[Mbhi]),str(EddR_array[EddRi]), str(ebvbbb_array[ebvi])] = np.log10(bbb_nu), bbb_Fnu_red        
        
        return BBBFdict_4plot, parameters_names

    elif modelsettings['BBB']=='D12_S':

        BBBFdict_4plot = dict()
        ## Call file containing all galaxy models     
        D12dict = cPickle.load(file(path + 'models/BBB/D12_S.pickle', 'rb'))    
        parameters_names =['logBHmass', 'logEddra']#D12dict['parameters']

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

        return BBBFdict_4plot, parameters_names

    elif modelsettings['BBB']=='D12_K':

        BBBFdict_4plot = dict()
        ## Call file containing all galaxy models     
        D12dict = cPickle.load(file(path + 'models/BBB/D12_K.pickle', 'rb'))    
        parameters_names =['logBHmass', 'logEddra']#D12dict['parameters']

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

        return BBBFdict_4plot, parameters_names

    else:
        print ' '
        print 'ERROR: The model with the name "'+modelsettings['BBB']+'" does not exist.'


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
            TORUSFdict_4plot[str(S04dict['Nh-values'][nhi])] = tor_nu0, renorm_template('TO',tor_Fnu0)

        ## Name the parameters that compose the keys of the dictionary: TORUSFdict_4plot[key]. 
        ## Add the names in the same order as their values are arranged in the dictionary key above.   
        parameters_names = ['Nh']

        return TORUSFdict_4plot, parameters_names

    elif modelsettings['TORUS']=='NK0': 
        
        TORUSFdict_4plot  = dict()

        NK0 = cPickle.load(file(path + 'models/TORUS/nenkova_v0.pickle', 'rb'))  #torus_object = cPickle.load(file(path + 'models/TORUS/silva_v1.pickle', 'rb')) 
        incl_idx=len(NK0.SED) #nhidx=len(torus_object.SED)
        #Construct dictionaries 
        for incl_i in range(incl_idx): #for nhi in range(nhidx):

            tor_nu0, tor_Fnu0 = NK0.wave[incl_i], NK0.SED[incl_i].squeeze() #tor_nu0, tor_Fnu0 = torus_object.wave[nhi], torus_object.SED[nhi].squeeze()
            TORUSFdict_4plot[str(NK0.incl[incl_i])] = tor_nu0, renorm_template('TO',tor_Fnu0) #TORUSFdict_4plot[str(torus_object.nh[nhi])] = tor_nu0, tor_Fnu0

        parameters_names = ['incl']

        return TORUSFdict_4plot, parameters_names   


def RADIO(LIR, conv_factor, sb_nu0, sb_Fnu0, radioexcess):

    q_IR_r14 = np.random.normal(2.64, 0.26,1)
    alpha_syn  = -0.75
    alpha_th = -0.1
    nth_frac =0.9

    L14 = 10**(np.log10(LIR*conv_factor)-q_IR_r14)/3.75e12#1.4e9 #to Wats
    nu_spacing= (np.log10(sb_nu0)[2]-np.log10(sb_nu0)[1])
    radio_points = (np.log10(sb_nu0)[0]-9)/nu_spacing
    radio_nu14= np.log10(1.4e9)

    radio_nu = np.arange(np.log10(sb_nu0)[0]- nu_spacing*int(radio_points),np.log10(sb_nu0)[0], nu_spacing)
    radio_nu2 = np.concatenate((radio_nu, np.log10(sb_nu0)[:np.argmax(sb_Fnu0)]))
    all_nu=  np.concatenate((radio_nu, np.log10(sb_nu0)))

    Lsb = np.concatenate((sb_Fnu0[0]*1e-4*np.ones(len(all_nu)-len(sb_Fnu0)),sb_Fnu0))

    Lsyn_0 = 10**(-1*alpha_syn* np.log10(1.4e9/10**radio_nu2[0])) * L14*(nth_frac)* rad_excess
    Lsyn_rad = Lsyn_0 * 10**(alpha_syn* np.log10(10**radio_nu2/10**radio_nu2[0])) 
    Lsyn = np.concatenate((Lsyn_rad, Lsyn_rad[-1]*1e-4*np.ones(len(all_nu)-len(Lsyn_rad))))

    Lff_0 = 10**(-1*alpha_th* np.log10(1.4e9/10**radio_nu2[0])) * L14*(1.-nth_frac)
    Lff_rad= Lff_0* 10**(alpha_th* np.log10(10**radio_nu2/10**radio_nu2[0])) 
    Lsyn = np.concatenate((Lsyn_rad, Lsyn_rad[-1]*1e-4*np.ones(len(all_nu)-len(Lsyn_rad))))

    Lir_rad= Lsb+Lsyn+Lff

    return  all_nu, Lir_rad



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
   
def fluxlambda_2_fluxnu (flux_lambda, wl_angst):

    c = 2.99792458e8

    flux_nu = flux_lambda * (wl_angst**2. ) / c /Angstrom
    return flux_nu

"""---------------------------------------------
             COMPUTED QUANTITIES
-----------------------------------------------"""

def stellar_info(chain, data):

    """
    computes stellar masses and SFRs
    """

    gal_obj,_,_,_ = data.dictkey_arrays
    _,_,_,_,_,_,_,_,SFRdict,_= data.dict_modelfluxes

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
            gal_obj.pick_nD(tuple([tau_mcmc[i], age_mcmc[i], 0.]))
            tau_dct, age_dct, ebvg_dct=gal_obj.matched_parkeys
            SFR_mcmc =SFRdict[tau_dct, age_dct]
        elif len(gal_obj.par_names)==4:
            gal_obj.pick_nD(tuple([metal_mcmc[i], tau_mcmc[i], age_mcmc[i], 0.]))
            metal_dct,tau_dct, age_dct, ebvg_dct=gal_obj.matched_parkeys
            SFR_mcmc =SFRdict[metal_dct,tau_dct, age_dct]

        # Calculate Mstar. BC03 templates are normalized to M* = 1 M_sun. 
        # Thanks to Kenneth Duncan, and his python version of BC03, smpy
        Mstar = np.log10(N * 1) 
        #Calculate SFR. output is in [Msun/yr]. 
        #print tau_mcmc[i], age_mcmc[i], (N * SFR_mcmc), np.log10(N), SFR_mcmc 
        SFR = N * SFR_mcmc
        SFR_list.append(SFR.value)    
        Mstar_list.append(Mstar)    

    return np.array(Mstar_list) , np.array(SFR_list)


def stellar_info_array(chain_flat, data, Nthin_compute):

    """
    computes arrays of stellar masses and SFRs
    """

    Ns, Npar = np.shape(chain_flat) 
    #chain_thinned = chain_flat[0:Ns:int(Ns/Nthin_compute),:]
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
        Fnu_norm = Fnu/1e60
        return Fnu_norm

