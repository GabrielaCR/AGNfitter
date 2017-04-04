
"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DICTIONARIES_AGNFitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains all functions which are needed to construct the total model of AGN. 

##For constructing a new dictionary, 
(in cases: 1)add a filter which is not included, 
2) need finer grid for better S/N data)
 see DICTIONARIES_AGNfitter.py  
    

"""
import numpy as np
import sys
from collections import defaultdict

import MODEL_AGNfitter as model
from scipy.integrate  import trapz
import time
import cPickle
import shelve
from astropy import units as u 


class MODELSDICT:


    """
    Class MODELSDICT

    Builds a dictionary of model templates. 

    ##input: 
    - filename of the dictionary you want to create
    - the path whre it will be located

    - Also variables self.ebvgal_array,self.ebvbbb_array, self.z_array
      can be change by the user, for a finer grid in this parameters.
    ##bugs: 

    """     

    def __init__(self, filename, path, filters):
        self.filename = filename
        self.path=path
        self.ebvgal_array = np.array(np.arange(0.,100.,5.)/100)
        self.ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)
        self.z_array = filters['dict_zarray']
        self.filterset = filters['Bandset']
        self.filters = filters

    def build(self):

        f = open(self.filename, 'wb')

        COSMOS_modelsdict = dict()

        print 'MODELSDICT.build'
        print 'Constructing Dictionary of models.' 
        print '--------------------------------------'
        print 'Make sure the filterset contains all the photometric bands'
        print 'needed by your catalog.'
        print 'This process might take a while, but you have to do it only once.'
        print 'If you interrupt it, please trash the empty file created.'
        print ''

        i=0
        dictionary_progressbar(i, len(self.z_array), prefix = 'Dict:', suffix = 'Complete', barLength = 50)
    
        for z in self.z_array:                
            i += 1
            filterdict = filter_dictionaries(self.filterset, self.path, self.filters)
            dict_modelsfiltered = self.construct_dictionaryarray_filtered(z, filterdict, self.path)
            COSMOS_modelsdict[str(z)] = dict_modelsfiltered

            time.sleep(0.01)
            dictionary_progressbar(i, len(self.z_array), prefix = 'Dict:', suffix = 'Complete', barLength = 50)



        print 'Dictionary has been created in :', self.filename

        cPickle.dump(COSMOS_modelsdict, f, protocol=2)
        f.close()



    def construct_dictionaryarray_filtered(self, z, filterdict,path):

        """
        Construct the dictionaries of fluxes at bands (to compare to data), 
        and dictionaries of fluxes over the whole spectrum, for plotting.
        """

        GALAXYFdict_filtered = dict()
        GALAXY_SFRdict = dict()
        STARBURSTFdict_filtered = dict()        
        BBBFdict_filtered = dict()
        TORUSFdict_filtered = dict()

        GALAXYFdict_4plot = dict()
        STARBURSTFdict_4plot = dict()        
        BBBFdict_4plot = dict()
        TORUSFdict_4plot = dict()



        #OPENING TEMPLATES AND BUILDING DICTIONARIES


        #Call object containing all galaxy models     
        galaxy_object = cPickle.load(file(path + 'models/GALAXY/bc03_275templates.pickle', 'rb')) 
        _, ageidx, tauidx, _, _,_ =  np.shape(galaxy_object.SED)
        #Construct dictionaries 
        for taui in range(tauidx):
            for agei in range(ageidx):
                t1= time.time()
                gal_wl, gal_Fwl =  galaxy_object.wave, galaxy_object.SED[:,agei,taui,:,:,:].squeeze()
                gal_nus= gal_wl.to(u.Hz, equivalencies=u.spectral())[::-1]#invert
                gal_Fnu= (gal_Fwl * 3.34e-19 * gal_wl**2.)[::-1]  
                gal_SFR= galaxy_object.SFR[:,agei,taui,:,:].squeeze()
                GALAXY_SFRdict[str(galaxy_object.tau.value[taui]),str(galaxy_object.tg.value[agei])] = gal_SFR

                for EBV_gal in self.ebvgal_array:
                    #Apply reddening            
                    gal_nu, gal_Fnu_red = model.GALAXYred_Calzetti(gal_nus.value[0:len(gal_nus):3], gal_Fnu.value[0:len(gal_nus):3], EBV_gal)                    
                    GALAXYFdict_4plot[str(galaxy_object.tau.value[taui]),str(galaxy_object.tg.value[agei]), str(EBV_gal)] = \
                                                                                            np.log10(gal_nu), gal_Fnu_red
                    #Projection of filter curves on models
                    bands,  gal_Fnu_filtered =  model.filters1(np.log10(gal_nu), gal_Fnu_red, filterdict, z)            
                    GALAXYFdict_filtered[str(galaxy_object.tau.value[taui]),str(galaxy_object.tg.value[agei]), str(EBV_gal)] = \
                                                                                                        bands, gal_Fnu_filtered
        #Call object containing all starburst models     
        starburst_object = cPickle.load(file(path + 'models/STARBURST/dalehelou_charyelbaz_v1.pickle', 'rb')) 
        irlumidx = len(starburst_object.SED)
        #Construct dictionaries 
        for irlumi in range(irlumidx):
            sb_nu0, sb_Fnu0 = starburst_object.wave[irlumi], starburst_object.SED[irlumi].squeeze()
            STARBURSTFdict_4plot[str(starburst_object.irlum[irlumi])] = sb_nu0, sb_Fnu0
            bands, sb_Fnu_filtered = model.filters1(sb_nu0, sb_Fnu0, filterdict, z)
            STARBURSTFdict_filtered[str(starburst_object.irlum[irlumi])] = bands, sb_Fnu_filtered
            if np.amax(sb_Fnu_filtered) == 0:
                print 'Error: something is wrong in the calculation of STARBURST flux'



        #No object to call since bbb is only one model     
        bbb_object = cPickle.load(file(path + 'models/BBB/richards.pickle', 'rb')) 

        bbb_nu, bbb_Fnu = bbb_object.wave, bbb_object.SED.squeeze()
        #Construct dictionaries
        for EBV_bbb in self.ebvbbb_array:
            bbb_nu0, bbb_Fnu_red = model.BBBred_Prevot(bbb_nu, bbb_Fnu, EBV_bbb, z )
            BBBFdict_4plot[str(EBV_bbb)] =bbb_nu0, bbb_Fnu_red
            bands, bbb_Fnu_filtered = model.filters1(bbb_nu0, bbb_Fnu_red, filterdict,z)
            BBBFdict_filtered[str(EBV_bbb)] = bands, bbb_Fnu_filtered
            if np.amax(bbb_Fnu_filtered) == 0:
                print 'Error: something is wrong in the calculation of BBB flux'            



        #Call object containing all torus models     
        torus_object = cPickle.load(file(path + 'models/TORUS/silva_v1.pickle', 'rb')) 
        nhidx=len(torus_object.SED)
        #Construct dictionaries 
        for nhi in range(nhidx):

            tor_nu0, tor_Fnu0 = torus_object.wave[nhi], torus_object.SED[nhi].squeeze()
            TORUSFdict_4plot[str(torus_object.nh[nhi])] = tor_nu0, tor_Fnu0

            bands, tor_Fnu_filtered = model.filters1(tor_nu0, tor_Fnu0, filterdict, z)
            TORUSFdict_filtered[str(torus_object.nh[nhi])] = bands, tor_Fnu_filtered
            if np.amax(tor_Fnu_filtered) == 0:
                print 'Error: something is wrong in the calculation of TORUS flux'




        return STARBURSTFdict_filtered , BBBFdict_filtered, GALAXYFdict_filtered, TORUSFdict_filtered, \
               STARBURSTFdict_4plot , BBBFdict_4plot, GALAXYFdict_4plot, TORUSFdict_4plot,GALAXY_SFRdict
               




def dictkey_arrays(MODELSdict):

    """
    Construct the dictionaries of fluxes at bands (to campare to data), 
    and dictionaries of fluxes over the whole spectrum, for plotting.

    ##input:

    ##output:
    """

    STARBURSTFdict , BBBFdict, GALAXYFdict, TORUSFdict, _,_,_,_,GALAXY_SFRdict= MODELSdict
    tau_dict= np.array(list(GALAXYFdict.keys()))[:,0]
    age_dict= np.array(list(GALAXYFdict.keys()))[:,1]
    ebvg_dict = np.array(list(GALAXYFdict.keys()))[:,2]

    irlum_dict = np.array(list(STARBURSTFdict.keys()))
    nh_dict = np.array(list(TORUSFdict.keys()))
    ebvb_dict = np.array(list(BBBFdict.keys()))


    #For computational reasons (to be used in PARAMETERspace_AGNfitter.py)
    class gal_class:
        def __init__(self, tau_dict, age_dict, ebvg_dict):
            self.tau_dict =tau_dict
            self.age_dict= age_dict
            self.ebvg_dict = ebvg_dict
            self.tau_dict_float =tau_dict.astype(float)
            self.age_dict_float= age_dict.astype(float)
            self.ebvg_dict_float = ebvg_dict.astype(float)

        def nearest_par2dict(self, tau, age, ebvg):    
            taui =np.abs(self.tau_dict_float-tau).argmin()
            agei= np.abs(self.age_dict_float-age).argmin()
            ebvgi = np.abs(self.ebvg_dict_float-ebvg).argmin()

            self.t = tau_dict[taui]
            self.a= age_dict[agei]
            self.e= ebvg_dict[ebvgi]

    gal_obj = gal_class(tau_dict, age_dict, ebvg_dict)

    return gal_obj, irlum_dict, nh_dict, ebvb_dict, GALAXY_SFRdict





def filter_dictionaries(filterset, path, filters):

    """
    Constructs the dictionaries of fluxes 
    1) specifically for your photometric bands (to campare to data), and
    2) dictionaries of fluxes for the whole spectrum, for plotting.

    input
    -------
    - filterset: Here we have two types of filterset: 
    'BANDSET_default' or 'BANDSET_settings'.
    This was specified from the RUN_AGNfitter_multi.py script.

    'BANDSET_default' includes bands needed for the example.
    'BANDSET_settings' includes all bands you specify in RUN_AGNfitter_multi.py. 

    dependency
    ----------
    This function is called in the CLASS MODELSDICT

    """
    H500band_file = path + 'models/FILTERS/HERSCHEL/SPIRE_500mu.txt'
    H500_lambda, H500_factor =  np.loadtxt(H500band_file, usecols=(0,1),unpack= True)

    H350band_file = path + 'models/FILTERS/HERSCHEL/SPIRE_350mu.txt'
    H350_lambda, H350_factor =  np.loadtxt(H350band_file, usecols=(0,1),unpack= True)

    H250band_file = path + 'models/FILTERS/HERSCHEL/SPIRE_250mu.txt'
    H250_lambda, H250_factor =  np.loadtxt(H250band_file, usecols=(0,1),unpack= True)

    H160band_file = path + 'models/FILTERS/HERSCHEL/PACS_160mu.txt'
    H160_lambda, H160_factor =  np.loadtxt(H160band_file, usecols=(0,1),unpack= True)

    H100band_file =path + 'models/FILTERS/HERSCHEL/PACS_100mu.txt'
    H100_lambda, H100_factor =  np.loadtxt(H100band_file, usecols=(0,1),unpack= True)

    #SPITZER
    M160band_file = path + 'models/FILTERS/SPITZER/mips160.res'
    M160_lambda, M160_factor =  np.loadtxt(M160band_file, usecols=(0,1),unpack= True)

    M70band_file = path + 'models/FILTERS/SPITZER/mips70.res'
    M70_lambda, M70_factor =  np.loadtxt(M70band_file, usecols=(0,1),unpack= True)

    M24band_file =  path + 'models/FILTERS/SPITZER/mips24.res'
    M24_lambda, M24_factor =  np.loadtxt(M24band_file, usecols=(0,1),unpack= True)

    #IRAC
    I4band_file = path + 'models/FILTERS/SPITZER/irac_ch4.res'
    I4_lambda, I4_factor =  np.loadtxt(I4band_file, usecols=(0,1),unpack= True)

    I3band_file =  path + 'models/FILTERS/SPITZER/irac_ch3.res'
    I3_lambda, I3_factor =  np.loadtxt(I3band_file, usecols=(0,1),unpack= True)

    I2band_file = path + 'models/FILTERS/SPITZER/irac_ch2.res'
    I2_lambda, I2_factor =  np.loadtxt(I2band_file, usecols=(0,1),unpack= True)

    I1band_file = path + 'models/FILTERS/SPITZER/irac_ch1.res'
    I1_lambda, I1_factor =  np.loadtxt(I1band_file, usecols=(0,1),unpack= True)

    #WISE
    W4band_file = path + 'models/FILTERS/WISE/NRSR-W4.txt'
    W4_lambda, W4_factor =  np.loadtxt(W4band_file, usecols=(0,1),unpack= True)

    W3band_file =  path + 'models/FILTERS/WISE/NRSR-W3.txt'
    W3_lambda, W3_factor =  np.loadtxt(W3band_file, usecols=(0,1),unpack= True)

    W2band_file = path + 'models/FILTERS/WISE/NRSR-W2.txt'
    W2_lambda, W2_factor =  np.loadtxt(W2band_file, usecols=(0,1),unpack= True)

    W1band_file = path + 'models/FILTERS/WISE/NRSR-W1.txt'
    W1_lambda, W1_factor =  np.loadtxt(W1band_file, usecols=(0,1),unpack= True)

    #2mass
    Kband_file = path + 'models/FILTERS/2MASS/Ks_2mass.res'
    K_lambda, K_factor =  np.loadtxt(Kband_file, usecols=(0,1),unpack= True)

    Hband_file = path + 'models/FILTERS/2MASS/H_2mass.res'
    H_lambda, H_factor =  np.loadtxt(Hband_file, usecols=(0,1),unpack= True)

    Jband_file = path + 'models/FILTERS/2MASS/J_2mass.res'
    J_lambda, J_factor =  np.loadtxt(Jband_file, usecols=(0,1),unpack= True)

    #VISTA

    Huvband_file = path + 'models/FILTERS/VISTA/H_uv.res'
    Huv_lambda, Huv_factor =  np.loadtxt(Huvband_file, usecols=(0,1),unpack= True)

    Juvband_file = path + 'models/FILTERS/VISTA/J_uv.res'
    Juv_lambda, Juv_factor =  np.loadtxt(Juvband_file, usecols=(0,1),unpack= True)

    Kuvband_file = path + 'models/FILTERS/VISTA/K_uv.res'
    Kuv_lambda, Kuv_factor =  np.loadtxt(Kuvband_file, usecols=(0,1),unpack= True)

    Yuvband_file = path + 'models/FILTERS/VISTA/Y_uv.res'
    Yuv_lambda, Yuv_factor =  np.loadtxt(Yuvband_file, usecols=(0,1),unpack= True)

    #CHFT ugriz
    uband_file_CHFT = path + 'models/FILTERS/CHFT/u_megaprime_sagem.res'
    u_lambda_CHFT, u_factor_CHFT =  np.loadtxt(uband_file_CHFT, usecols=(0,1),unpack= True)

    gband_file_CHFT = path + 'models/FILTERS/CHFT/g_megaprime_sagem.res'
    g_lambda_CHFT, g_factor_CHFT =  np.loadtxt(gband_file_CHFT, usecols=(0,1),unpack= True)

    rband_file_CHFT = path + 'models/FILTERS/CHFT/r_megaprime_sagem.res'
    r_lambda_CHFT, r_factor_CHFT =  np.loadtxt(rband_file_CHFT, usecols=(0,1),unpack= True)

    iband_file_CHFT = path + 'models/FILTERS/CHFT/i_megaprime_sagem.res'
    i_lambda_CHFT, i_factor_CHFT =  np.loadtxt(iband_file_CHFT, usecols=(0,1),unpack= True)

    zband_file_CHFT = path + 'models/FILTERS/CHFT/z_megaprime_sagem.res'
    z_lambda_CHFT, z_factor_CHFT =  np.loadtxt(zband_file_CHFT, usecols=(0,1),unpack= True)

    #SDSS ugriz
    uband_file_SDSS = path + 'models/FILTERS/SDSS/u_SDSS.res'
    u_lambda_SDSS, u_factor_SDSS =  np.loadtxt(uband_file_SDSS, usecols=(0,1),unpack= True)

    gband_file_SDSS = path + 'models/FILTERS/SDSS/g_SDSS.res'
    g_lambda_SDSS, g_factor_SDSS =  np.loadtxt(gband_file_SDSS, usecols=(0,1),unpack= True)

    rband_file_SDSS = path + 'models/FILTERS/SDSS/r_SDSS.res'
    r_lambda_SDSS, r_factor_SDSS =  np.loadtxt(rband_file_SDSS, usecols=(0,1),unpack= True)

    iband_file_SDSS = path + 'models/FILTERS/SDSS/i_SDSS.res'
    i_lambda_SDSS, i_factor_SDSS =  np.loadtxt(iband_file_SDSS, usecols=(0,1),unpack= True)

    zband_file_SDSS = path + 'models/FILTERS/SDSS/z_SDSS.res'
    z_lambda_SDSS, z_factor_SDSS =  np.loadtxt(zband_file_SDSS, usecols=(0,1),unpack= True)

    #SUBARU

    gband_file =path + 'models/FILTERS/SUBARU/g_subaru.res'
    g_lambda,g_factor =  np.loadtxt(gband_file, usecols=(0,1),unpack= True)

    rband_file = path + 'models/FILTERS/SUBARU/r_subaru.res'
    r_lambda,r_factor =  np.loadtxt(rband_file, usecols=(0,1),unpack= True)

    iband_file = path + 'models/FILTERS/SUBARU/i_subaru.res'
    i_lambda,i_factor =  np.loadtxt(iband_file, usecols=(0,1),unpack= True)

    zband_file =path + 'models/FILTERS/SUBARU/z_subaru.res'
    z_lambda, z_factor =  np.loadtxt(zband_file, usecols=(0,1),unpack= True)

    Bband_file = path + 'models/FILTERS/SUBARU/B_subaru.res'
    B_lambda, B_factor =  np.loadtxt(Bband_file, usecols=(0,1),unpack= True)

    Vband_file = path + 'models/FILTERS/SUBARU/V_subaru.res'
    V_lambda, V_factor =  np.loadtxt(Vband_file, usecols=(0,1),unpack= True)

    #GALEX
    NUVband_file = path + 'models/FILTERS/GALEX/galex2500.res'
    NUV_lambda, NUV_factor =  np.loadtxt(NUVband_file, usecols=(0,1),unpack= True)

    FUVband_file = path + 'models/FILTERS/GALEX/galex1500.res'
    FUV_lambda, FUV_factor =  np.loadtxt(FUVband_file, usecols=(0,1),unpack= True)


                
    if filterset == 'BANDSET_default':


        #List of file names
        files = [ H500band_file, H350band_file, H250band_file, M24band_file, I4band_file ,\
         I3band_file, I2band_file, I1band_file, Kband_file, Hband_file, Jband_file,\
          Yuvband_file, zband_file , iband_file, rband_file,  Bband_file,  uband_file_CHFT, \
          NUVband_file]

        #List of all lambdas
        lambdas = [H500_lambda, H350_lambda, H250_lambda, M24_lambda, I4_lambda , I3_lambda, \
        I2_lambda, I1_lambda,  K_lambda, H_lambda, J_lambda, Yuv_lambda,  z_lambda, i_lambda, \
        r_lambda, B_lambda,  u_lambda_CHFT, NUV_lambda]

        #list of all factors corresponding to the lambdas
        factors = [ H500_factor, H350_factor, H250_factor, M24_factor, I4_factor , I3_factor, \
        I2_factor, I1_factor, K_factor, H_factor, J_factor, Yuv_factor, z_factor, i_factor, \
        r_factor,  B_factor,  u_factor_CHFT, NUV_factor]


        #dictionaries lambdas_dict, factors_dict
        files_dict = defaultdict(list)
        lambdas_dict = defaultdict(list)
        factors_dict = defaultdict(list)
        central_nu_list=[]

        for i in range(len(files)):

            c=    2.997e8
            Angstrom = 1e10

            central_lamb = np.sum(lambdas[i]*factors[i])/np.sum(factors[i])
            central_nu = float(np.log10((Angstrom*c)/central_lamb))

            files_dict[central_nu].append(files[i])
            lambdas_dict[central_nu].append(lambdas[i])
            factors_dict[central_nu].append(factors[i])

            central_nu_list.append(central_nu)


    if filterset == 'BANDSET_settings':

        files =[]
        lambdas = []
        factors = []

        if filters['SPIRE500']:
            files.append(H500band_file)
            lambdas.append(H500_lambda)
            factors.append(H500_factor)
        if filters['SPIRE350']:
            files.append(H350band_file)
            lambdas.append(H350_lambda)
            factors.append(H350_factor)
        if filters['SPIRE250']:
            files.append(H250band_file)
            lambdas.append(H250_lambda)
            factors.append(H250_factor)
        if filters['PACS160']:
            files.append(H160band_file)
            lambdas.append(H160_lambda)
            factors.append(H160_factor)
        if filters['PACS100']:
            files.append(H100band_file)
            lambdas.append(H100_lambda)
            factors.append(H100_factor)

        if filters['MIPS160']:      
            files.append(M160band_file)
            lambdas.append(M160_lambda)
            factors.append(M160_factor)
        if filters['MIPS70']:    
            files.append(M70band_file)
            lambdas.append(M70_lambda)
            factors.append(M70_factor)
        if filters['MIPS24']:
            files.append(M24band_file)
            lambdas.append(M24_lambda)
            factors.append(M24_factor)

        if filters['IRAC4']:   
            files.append(I4band_file)
            lambdas.append(I4_lambda)
            factors.append(I4_factor)         
        if filters['IRAC3']:
            files.append(I3band_file)
            lambdas.append(I3_lambda)
            factors.append(I3_factor)    
        if filters['IRAC2']:
            files.append(I2band_file)
            lambdas.append(I2_lambda)
            factors.append(I2_factor)                
        if filters['IRAC1']:
            files.append(I1band_file)
            lambdas.append(I1_lambda)
            factors.append(I1_factor)    

        if filters['WISE4']:
            files.append(W4band_file)
            lambdas.append(W4_lambda)
            factors.append(W4_factor)    
        if filters['WISE3']:
            files.append(W3band_file)
            lambdas.append(W3_lambda)
            factors.append(W3_factor)    
        if filters['WISE2']:
            files.append(W2band_file)
            lambdas.append(W2_lambda)
            factors.append(W2_factor)    
        if filters['WISE1']:
            files.append(W1band_file)
            lambdas.append(W1_lambda)
            factors.append(W1_factor)    

        if filters['Ks_2mass']:
            files.append(Kband_file)
            lambdas.append(K_lambda)
            factors.append(K_factor)    
        if filters['H_2mass']:
            files.append(Hband_file)
            lambdas.append(H_lambda)
            factors.append(H_factor)    
        if filters['J_2mass']:
            files.append(Jband_file)
            lambdas.append(J_lambda)
            factors.append(J_factor)    

        if filters['H_VISTA']:
            files.append(Huvband_file)
            lambdas.append(Huv_lambda)
            factors.append(Huv_factor)    
        if filters['J_VISTA']:
            files.append(Juvband_file)
            lambdas.append(Juv_lambda)
            factors.append(Juv_factor)    
        if filters['K_VISTA']:
            files.append(Kuvband_file)
            lambdas.append(Kuv_lambda)
            factors.append(Kuv_factor)        
        if filters['Y_VISTA']:
            files.append(Yuvband_file)
            lambdas.append(Yuv_lambda)
            factors.append(Yuv_factor)    

        if filters['g_SUBARU']:
            files.append(gband_file)
            lambdas.append(g_lambda)
            factors.append(g_factor)    
        if filters['r_SUBARU']:
            files.append(rband_file)
            lambdas.append(r_lambda)
            factors.append(r_factor)    
        if filters['i_SUBARU']:  
            files.append(iband_file)
            lambdas.append(i_lambda)
            factors.append(i_factor)    
        if filters['z_SUBARU']:
            files.append(zband_file)
            lambdas.append(z_lambda)
            factors.append(z_factor)    
        if filters['B_SUBARU']:
            files.append(Bband_file)
            lambdas.append(B_lambda)
            factors.append(B_factor)    
        if filters['V_SUBARU']:  
            files.append(Vband_file)
            lambdas.append(V_lambda)
            factors.append(V_factor)    

        if filters['u_CHFT']:  
            files.append(uband_file_CHFT)
            lambdas.append(u_lambda_CHFT)
            factors.append(u_factor_CHFT)    
        if filters['g_CHFT']:
            files.append(gband_file_CHFT)
            lambdas.append(g_lambda_CHFT)
            factors.append(g_factor_CHFT)    
        if filters['r_CHFT']:
            files.append(rband_file_CHFT)
            lambdas.append(r_lambda_CHFT)
            factors.append(r_factor_CHFT)    
        if filters['i_CHFT']:  
            files.append(iband_file_CHFT)
            lambdas.append(i_lambda_CHFT)
            factors.append(i_factor_CHFT)    
        if filters['z_CHFT']:
            files.append(zband_file_CHFT)
            lambdas.append(z_lambda_CHFT)
            factors.append(z_factor_CHFT)    

        if filters['u_SDSS']:  
            files.append(uband_file_SDSS)
            lambdas.append(u_lambda_SDSS)
            factors.append(u_factor_SDSS)    
        if filters['g_SDSS']:
            files.append(gband_file_SDSS)
            lambdas.append(g_lambda_SDSS)
            factors.append(g_factor_SDSS)    
        if filters['r_SDSS']:
            files.append(rband_file_SDSS)
            lambdas.append(r_lambda_SDSS)
            factors.append(r_factor_SDSS)    
        if filters['i_SDSS']:  
            files.append(iband_file_SDSS)
            lambdas.append(i_lambda_SDSS)
            factors.append(i_factor_SDSS)    
        if filters['z_SDSS']:
            files.append(zband_file_SDSS)
            lambdas.append(z_lambda_SDSS)
            factors.append(z_factor_SDSS)    

        if filters['GALEX_2500']:
            files.append(NUVband_file)
            lambdas.append(NUV_lambda)
            factors.append(NUV_factor)    
        if filters['GALEX_1500']:
            files.append(FUVband_file)
            lambdas.append(FUV_lambda)
            factors.append(FUV_factor)    


        # make dictionaries lambdas_dict, factors_dict
        files_dict = defaultdict(list)
        lambdas_dict = defaultdict(list)
        factors_dict = defaultdict(list)
        central_nu_list=[]        

        for i in range(len(files)):

            c=    2.997e8
            Angstrom = 1e10

            central_lamb = np.sum(lambdas[i]*factors[i])/np.sum(factors[i])
            central_nu = float(np.log10((Angstrom*c)/central_lamb))

            files_dict[central_nu].append(files[i])
            lambdas_dict[central_nu].append(lambdas[i])
            factors_dict[central_nu].append(factors[i])

            central_nu_list.append(central_nu)

        central_nu_list = sorted(central_nu_list)
            
    return central_nu_list, files_dict, lambdas_dict, factors_dict



def dictionary_progressbar (iteration, total, prefix = '', suffix = '', decimals = 2, barLength = 100):
    """
    Print progress bar of dictionary construction
    """
    filledLength    = int(round(barLength * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '>' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()


