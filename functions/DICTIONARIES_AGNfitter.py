
"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DICTIONARIES_AGNFitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains all functions which are needed to construct the total model of AGN. 

##For constructing a new dictionary, 
(in cases: 1)add a filter which is not included, 
2) need finer grid for better S/N data)
 see DICTIONARIES_AGNfitter.py  
    

"""
import sys,os
import numpy as np
import sys
from collections import defaultdict

import MODEL_AGNfitter as model
import FILTERS_AGNfitter as filterpy

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
    - the AGNfitter path is in your computer 
    - the filters settings (dictionary from the settings file)

    - Also variables self.ebvgal_array,self.ebvbbb_array, self.z_array
      can be changed by the user, for a finer grid in this parameters.
    ##bugs: 

    """     

    def __init__(self, filename, path, filters):
        self.filename = filename
        self.path=path
        self.ebvgal_array = np.array(np.arange(0.,100.,5.)/100)
        self.ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)

        ## To be called form filters
        self.z_array = filters['dict_zarray']
        self.filters_list = filters
        if os.path.lexists(filename):
            self.fo = cPickle.load(file(filename, 'rb'))
            self.filters = self.fo.filternames
            self.filterset_name = self.fo.name
        else:
            self.fo = filterpy.create_filtersets(filters, path)
            self.filters = self.fo.filternames
            self.filterset_name = self.fo.name

    def build(self):

        Modelsdict = dict()

        i=0
        dictionary_progressbar(i, len(self.z_array), prefix = 'Dict:', suffix = 'Complete', barLength = 50)
    
        for z in self.z_array:                
            i += 1
            #filterdict = filter_dictionaries(self.filterset, self.path, self.filters_list)
            filterdict = [self.fo.central_nu_array, self.fo.lambdas_dict, self.fo.factors_dict]
            dict_modelsfiltered = self.construct_dictionaryarray_filtered(z, filterdict, self.path)
            Modelsdict[str(z)] = dict_modelsfiltered
            time.sleep(0.01)
            dictionary_progressbar(i, len(self.z_array), prefix = 'Dict:', suffix = 'Complete', barLength = 50)

        self.MD = Modelsdict




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
    Summarizes the model dictionary keys.
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


