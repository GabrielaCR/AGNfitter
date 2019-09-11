
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
from scipy.interpolate import interp1d
import time
import cPickle
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

    def __init__(self, filename, path, filters, models):
        self.filename = filename
        self.path=path
        self.modelsettings=models

        ## To be called form filters
        self.z_array = filters['dict_zarray']
                    
        a = dict.fromkeys(filters)

        for i in range(len(a.keys())):

            if a.keys()[i] == 'add_filters' or a.keys()[i] == 'dict_zarray' or a.keys()[i] == 'add_filters_dict' or a.keys()[i] == 'path':

                a[a.keys()[i]] = filters[a.keys()[i]]               
                
            else:
                a[a.keys()[i]] = filters[a.keys()[i]]#[0]

        #print a
        self.filters_list = a
        if os.path.lexists(filename):
            self.fo = cPickle.load(file(filename, 'rb'))
            self.filters = self.fo.filternames
            self.filterset_name = self.fo.name
        else:
            self.fo = filterpy.create_filtersets(a, path)
            self.filters = self.fo.filternames
            self.filterset_name = self.fo.name

    def build(self):

        Modelsdict = dict()

        i=0
        dictionary_progressbar(i, len(self.z_array), prefix = 'Dict:', suffix = 'Complete', barLength = 50)
    
        for z in self.z_array:                
            i += 1
            filterdict = [self.fo.central_nu_array, self.fo.lambdas_dict, self.fo.factors_dict]
            dict_modelsfiltered = construct_dictionaryarray_filtered(z, filterdict, self.path, self.modelsettings)
            Modelsdict[str(z)] = dict_modelsfiltered
            time.sleep(0.01)
            dictionary_progressbar(i, len(self.z_array), prefix = 'Dict:', suffix = 'Complete', barLength = 50)

        self.MD = Modelsdict


def construct_dictionaryarray_filtered( z, filterdict,path, modelsettings):

    """
    Construct the dictionaries of fluxes at bands (to compare to data), 
    and dictionaries of fluxes over the whole spectrum, for plotting.
    All calculations are done at one given redshift.
    """

    GALAXYFdict_filtered = dict()
    STARBURSTFdict_filtered = dict()        
    BBBFdict_filtered = dict()
    TORUSFdict_filtered = dict()

  
    GALAXYFdict_4plot, GALAXY_SFRdict, galaxy_parnames  = model.GALAXY(path, modelsettings)
    for c in GALAXYFdict_4plot.keys():
                gal_nu, gal_Fnu=GALAXYFdict_4plot[c]               
                bands,  gal_Fnu_filtered =  filtering_models(gal_nu, gal_Fnu, filterdict, z)            
                GALAXYFdict_filtered[c] = bands, gal_Fnu_filtered

    STARBURSTFdict_4plot, starburst_parnames  = model.STARBURST(path, modelsettings)
    for c in STARBURSTFdict_4plot.keys():
                sb_nu, sb_Fnu=STARBURSTFdict_4plot[c]               
                bands, sb_Fnu_filtered  =  filtering_models(sb_nu, sb_Fnu, filterdict, z)            
                STARBURSTFdict_filtered[c] = bands, sb_Fnu_filtered

    BBBFdict_4plot, bbb_parnames = model.BBB(path, modelsettings)
    for c in BBBFdict_4plot.keys():
                bbb_nu, bbb_Fnu=BBBFdict_4plot[c]               
                bands,  bbb_Fnu_filtered =  filtering_models(bbb_nu, bbb_Fnu, filterdict, z)            
                BBBFdict_filtered[c] = bands, bbb_Fnu_filtered

    TORUSFdict_4plot, torus_parnames  = model.TORUS(path, modelsettings)
    for c in TORUSFdict_4plot.keys():
                tor_nu, tor_Fnu=TORUSFdict_4plot[c]               
                bands, tor_Fnu_filtered  =  filtering_models(tor_nu, tor_Fnu, filterdict, z)            
                TORUSFdict_filtered[c] = bands, tor_Fnu_filtered

    norm_parnames = ['GA', 'SB', 'BB', 'TO' ]
    all_parnames = [galaxy_parnames, starburst_parnames,torus_parnames, bbb_parnames, norm_parnames]

    return STARBURSTFdict_filtered , BBBFdict_filtered, GALAXYFdict_filtered, TORUSFdict_filtered, \
           STARBURSTFdict_4plot , BBBFdict_4plot, GALAXYFdict_4plot, TORUSFdict_4plot,GALAXY_SFRdict, all_parnames
           

def dictkey_arrays(MODELSdict):

    """
    Summarizes the model dictionary keys and does the interpolation to nearest value in grid.
    used to be transporte to data

    ##input:

    ##output:
    """

    STARBURSTFdict , BBBFdict, GALAXYFdict, TORUSFdict, _,_,_,_,GALAXY_SFRdict, all_parnames= MODELSdict

    galaxy_parkeys= np.array(list(GALAXYFdict.keys()))
    starburst_parkeys = np.array(list(STARBURSTFdict.keys()))
    torus_parkeys = np.array(list(TORUSFdict.keys()))
    bbb_parkeys = np.array(list(BBBFdict.keys()))

    class pick_obj:
            def __init__(self, par_names,pars_modelkeys):

                self.pars_modelkeys=pars_modelkeys.T
                self.pars_modelkeys_float =self.pars_modelkeys.astype(float)
                self.par_names = par_names

            def pick_nD(self, pars_mcmc): 
                self.matched_parkeys = []
                if len(pars_mcmc)==1:
                    for i in range(len(pars_mcmc)):   
                        matched_idx =np.abs(self.pars_modelkeys_float-pars_mcmc[i]).argmin()
                        matched_parkey = self.pars_modelkeys[matched_idx]
                        self.matched_parkeys =matched_parkey
                else:
                    for i in range(len(pars_mcmc)):   
                        matched_idx =np.abs(self.pars_modelkeys_float[i]-pars_mcmc[i]).argmin()
                        matched_parkey = self.pars_modelkeys[i][matched_idx]
                        self.matched_parkeys.append(matched_parkey)
                    self.matched_parkeys=tuple(self.matched_parkeys)

            def pick_1D(self, *pars_mcmc): 
                matched_idx =np.abs(self.pars_modelkeys_float-pars_mcmc).argmin()
                self.matched_parkeys = self.pars_modelkeys[matched_idx]

    galaxy_parnames, starburst_parnames,torus_parnames, bbb_parnames, norm_parnames = all_parnames
    gal_obj =pick_obj(galaxy_parnames,galaxy_parkeys)
    sb_obj =pick_obj(starburst_parnames,starburst_parkeys)
    tor_obj=pick_obj(torus_parnames,torus_parkeys)
    bbb_obj=pick_obj(bbb_parnames,bbb_parkeys)

    return gal_obj,sb_obj,tor_obj, bbb_obj 




def filtering_models( model_nus, model_fluxes, filterdict, z ):    

    """
    Projects the model SEDs into the filter curves of each photometric band.

    ##input:
    - model_nus: template frequencies [log10(nu)]
    - model_fluxes: template fluxes [F_nu]
    - filterdict: dictionary with all band filter curves' information.
                  To change this, add one band and filter curve, etc,
                  look at DICTIONARIES_AGNfitter.py
    - z: redshift

    ##output:
    - bands [log10(nu)]
    - Filtered fluxes at these bands [F_nu]
    """

    bands, lambdas_dict, factors_dict = filterdict
    filtered_model_Fnus = []


    # Costumize model frequencies and fluxes [F_nu]
    # to same units as filter curves (to wavelengths [angstrom] and F_lambda)
    model_lambdas = nu2lambda_angstrom(model_nus) * (1+z)
    model_lambdas =  model_lambdas[::-1]
    model_fluxes_nu =  model_fluxes[::-1]
    model_fluxes_lambda = fluxnu_2_fluxlambda(model_fluxes_nu, model_lambdas) 
    mod2filter_interpol = interp1d(model_lambdas, model_fluxes_lambda, kind = 'nearest', bounds_error=False, fill_value=0.)            

    # For filter curve at each band. 
    # (Vectorised integration was not possible -> different filter-curve-arrays' sizes)
    for iband in bands:

        # Read filter curves info for each data point 
        # (wavelengths [angstrom] and factors [non])
        lambdas_filter = np.array(lambdas_dict[iband])
        factors_filter = np.array(factors_dict[iband])
        iband_angst = nu2lambda_angstrom(iband)

        # Interpolate the model fluxes to 
        #the exact wavelengths of filter curves
        modelfluxes_at_filterlambdas = mod2filter_interpol(lambdas_filter)
        # Compute the flux ratios, equivalent to the filtered fluxes: 
        # F = int(model)/int(filter)
        integral_model = trapz(modelfluxes_at_filterlambdas*factors_filter, x= lambdas_filter)
        integral_filter = trapz(factors_filter, x= lambdas_filter)     
        filtered_modelF_lambda = (integral_model/integral_filter)

        # Convert all from lambda, F_lambda  to Fnu and nu    
        filtered_modelFnu_atfilter_i = fluxlambda_2_fluxnu(filtered_modelF_lambda, iband_angst)
        filtered_model_Fnus.append(filtered_modelFnu_atfilter_i)

    return bands, np.array(filtered_model_Fnus)



## ---------------------------------------------------
c = 2.997e8
Angstrom = 1.e10

def fluxlambda_2_fluxnu (flux_lambda, wl_angst):
    """
    Calculate F_nu from F_lambda.
    """
    flux_nu = flux_lambda * (wl_angst**2. ) / c /Angstrom
    return flux_nu
def fluxnu_2_fluxlambda (flux_nu, wl_angst):
    """
    Calculate F_lambda from  F_nu.
    """
    flux_lambda = flux_nu / wl_angst**2 *c * Angstrom
    return flux_lambda #in angstrom
def nu2lambda_angstrom(nus):
    """
    Calculate wavelength [angstrom] from frequency [log Hz].
    """
    lambdas = c / (10**nus) * Angstrom
    return lambdas
## ------------------------------------------------





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


