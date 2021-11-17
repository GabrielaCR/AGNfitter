
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
from . import MODEL_AGNfitter as model
from . import FILTERS_AGNfitter as filterpy
from scipy.integrate  import trapz
from scipy.interpolate import interp1d
import time
import pickle 


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

    def __init__(self, filename, path, filters, models, nRADdata, nXRaysdata):
        self.filename = filename
        self.path=path
        self.modelsettings=models
        self.nRADdata = nRADdata  #Number of valid radio data
        self.nXRaysdata = nXRaysdata  #Number of valid Xrays data

        ## To be called form filters
        self.z_array = filters['dict_zarray']
                    
        a = dict.fromkeys(filters)

        for i in range(len(list(a.keys()))):

            if list(a.keys())[i] == 'add_filters' or list(a.keys())[i] == 'dict_zarray' or list(a.keys())[i] == 'add_filters_dict' or list(a.keys())[i] == 'path':

                a[list(a.keys())[i]] = filters[list(a.keys())[i]]               
                
            else:
                a[list(a.keys())[i]] = filters[list(a.keys())[i]]

        self.filters_list = a
        if os.path.lexists(filename):
            self.fo = pickle.load(open(filename, 'rb'))
            self.filters = self.fo.filternames
            self.filterset_name = self.fo.name
        else:
            self.fo = filterpy.create_filtersets(a, path)
            self.filters = self.fo.filternames
            self.filterset_name = self.fo.name

    def construct_dictionaryarray_filtered(self, z, filterdict):

        """
        Construct the dictionaries of fluxes at bands (to compare to data), 
        and dictionaries of fluxes over the whole spectrum, for plotting.
        All calculations are done at one given redshift.
        """

        self.GALAXYFdict = dict()
        self.STARBURSTFdict = dict()        
        self.BBBFdict = dict()
        self.TORUSFdict = dict()
        self.z= z
        self.filterdict= filterdict
        self.GALAXYFdict_4plot, self.GALAXY_SFRdict, self.GALAXYatt_dict, galaxy_parnames, galaxy_partypes,self.GALAXYfunctions  = model.GALAXY(self.path, self.modelsettings)
        for c in self.GALAXYFdict_4plot.keys():
                    gal_nu, gal_Fnu=self.GALAXYFdict_4plot[c]               
                    bands,  gal_Fnu_filtered =  filtering_models(gal_nu, gal_Fnu, filterdict, z)            
                    self.GALAXYFdict[c] = bands, gal_Fnu_filtered.flatten()

        self.STARBURSTFdict_4plot, self.STARBURST_LIRdict, starburst_parnames, starburst_partypes ,self.STARBURSTfunctions = model.STARBURST(self.path, self.modelsettings)
        for c in self.STARBURSTFdict_4plot.keys():
                    sb_nu, sb_Fnu=self.STARBURSTFdict_4plot[c]               
                    bands, sb_Fnu_filtered  =  filtering_models(sb_nu, sb_Fnu, filterdict, z)            
                    self.STARBURSTFdict[c] = bands, sb_Fnu_filtered.flatten()

        self.BBBFdict_4plot, bbb_parnames, bbb_partypes ,self.BBBfunctions= model.BBB(self.path, self.modelsettings, self.nXRaysdata)
        for c in self.BBBFdict_4plot.keys():
                    bbb_nu, bbb_Fnu=self.BBBFdict_4plot[c]             
                    bands,  bbb_Fnu_filtered =  filtering_models(bbb_nu, bbb_Fnu, filterdict, z)              
                    self.BBBFdict[c] = bands, bbb_Fnu_filtered.flatten()

        self.TORUSFdict_4plot, torus_parnames, torus_partypes ,self.TORUSfunctions = model.TORUS(self.path, self.modelsettings)
        for c in self.TORUSFdict_4plot.keys():
                    tor_nu, tor_Fnu=self.TORUSFdict_4plot[c]               
                    bands, tor_Fnu_filtered  =  filtering_models(tor_nu, tor_Fnu, filterdict, z)            
                    self.TORUSFdict[c] = bands, tor_Fnu_filtered.flatten()

        norm_parnames = ['GA', 'SB', 'BB', 'TO' ]
        norm_partypes = ['free', 'free', 'free', 'free' ]
        self.all_parnames = [galaxy_parnames, starburst_parnames,torus_parnames, bbb_parnames]
        self.all_partypes = [galaxy_partypes, starburst_partypes,torus_partypes, bbb_partypes]


        if self.modelsettings['RADIO']== True:  #If there are radio data available, the SEDfitting consider 5 components (AGN radio is the 5th)
            self.AGN_RADFdict = dict()
            self.AGN_RADFdict_4plot, agnrad_parnames, agnrad_partypes, self.AGN_RADfunctions  = model.AGN_RAD(self.path, self.modelsettings, self.nRADdata)
            for c in self.AGN_RADFdict_4plot.keys():
                agnrad_nu, agnrad_Fnu = self.AGN_RADFdict_4plot[c]               
                bands, agnrad_Fnu_filtered  =  filtering_models(agnrad_nu, agnrad_Fnu, filterdict, z)            
                self.AGN_RADFdict[c] = bands, agnrad_Fnu_filtered.flatten()
            norm_parnames.append('RAD')
            norm_partypes.append('free')
            self.all_parnames.extend((agnrad_parnames, norm_parnames))
            self.all_partypes.extend((agnrad_partypes, norm_partypes))

        else:
            self.all_parnames.append(norm_parnames)
            self.all_partypes.append(norm_partypes)


    def build(self):

        Modelsdict = dict()

        i=0
        dictionary_progressbar(i, len(self.z_array), prefix = 'Dict:', suffix = 'Complete', barLength = 50)
    
        for z in self.z_array:                
            i += 1
            filterdict = [self.fo.central_nu_array, self.fo.lambdas_dict, self.fo.factors_dict]
            dict_modelsfiltered = self.construct_dictionaryarray_filtered(z, filterdict)          #Models SEDs are redshifted and filtered
            Modelsdict[str(z)] = dict_modelsfiltered
            time.sleep(0.01)
            dictionary_progressbar(i, len(self.z_array), prefix = 'Dict:', suffix = 'Complete', barLength = 50)

        self.MD = Modelsdict

def dictkey_arrays(MD):

    """
    Summarizes the model dictionary keys and does the interpolation to nearest value in grid.
    used to be transporte to data

    ##input:

    ##output:
    """

    galaxy_parkeys= np.array(list(MD.GALAXYFdict.keys()))
    starburst_parkeys = np.array(list(MD.STARBURSTFdict.keys()))
    torus_parkeys = np.array(list(MD.TORUSFdict.keys()))
    bbb_parkeys = np.array(list(MD.BBBFdict.keys()))

    class get_model:
            def __init__(self, par_names, par_types, pars_modelkeys, modelsdict, z, functionidxs, functions):

                self.pars_modelkeys=pars_modelkeys.T
                self.pars_modelkeys_float =self.pars_modelkeys.astype(float)
                self.par_names = par_names
                self.par_types = par_types
                self.modelsdict = modelsdict
                self.functions = functions
                self.functionidxs=functionidxs
                self.z= z

            def pick_nD(self, pars_mcmc): 
                self.matched_parkeys = []
                self.matched_parkeys_grid = []

                if len(pars_mcmc)==1:
                    for i in range(len(pars_mcmc)):   
                        if self.par_types[i] == 'grid':
                            matched_idx =np.abs(self.pars_modelkeys_float-pars_mcmc[i]).argmin()  #Choose the parameter value closest to that found by mcmc
                            matched_parkey = self.pars_modelkeys[matched_idx]
                            self.matched_parkeys =matched_parkey
                            self.matched_parkeys_grid=self.matched_parkeys
                        elif self.par_types[i] == 'free':
                            self.matched_parkeys=pars_mcmc[i]                                    #Values found by mcmc in case of free parameters
                            self.matched_parkeys_grid=self.pars_modelkeys[0]
                        else: 
                            print('Error DICTIONARIES_AGNfitter.py: parameter type ',self.par_types, ' is unknown.')

                else:
                    for i in range(len(pars_mcmc)):
                        if self.par_types[i] == 'grid': 
                            matched_idx =np.abs(self.pars_modelkeys_float[i]-pars_mcmc[i]).argmin() #Choose the parameter value closest to that found by mcmc
                            matched_parkey = self.pars_modelkeys[i][matched_idx]
                            self.matched_parkeys.append(matched_parkey)
                            self.matched_parkeys_grid.append(matched_parkey)     
                        elif self.par_types[i] == 'free':
                            self.matched_parkeys.append(pars_mcmc[i])                               #Values found by mcmc in case of free parameters
                            self.matched_parkeys_grid.append(self.pars_modelkeys[i,0]) 
                        else: 
                            print('Error DICTIONARIES_AGNfitter.py: parameter type ',self.par_types, ' is unknown.')

                    self.matched_parkeys=tuple(self.matched_parkeys)

            def get_fluxes(self,  matched_parkeys):  #From the dictionary of redshifted and filtered models
                    
                    if 'free' not in self.par_types:   
                        return self.modelsdict[matched_parkeys]


                    elif 'free' in self.par_types and self.par_names[-1] == 'EBVgal':  #This is for the case of EBVgal == free
                        fcts=self.functions()
                        idxs=0
                        f=fcts[self.functionidxs[idxs]]
                        bands, Fnu = self.modelsdict[tuple(self.matched_parkeys_grid)] 
                        rest_bands = bands + np.log10((1+self.z))                        #Rest frame frequency
                        bandsf, Fnuf = f(10**rest_bands, Fnu, matched_parkeys[-1])       #Calzetti function need normal frequency (not log)
                        bandsf = np.log10(bandsf) - np.log10((1+self.z))                 #Come back to frequency corrected by redshift
                        return bandsf, Fnuf

                    elif 'free' in self.par_types and self.par_names[-1] == 'EBVbbb':  #This is for the case of EBVbbb == free
                        fcts=self.functions()
                        idxs=0
                        f=fcts[self.functionidxs[idxs]]
                        #R06 without X-Rays only have 1 parameter (EBV_bbb) and not a list of parameters so tuple() produce problems
                        if type(self.matched_parkeys_grid) != list:         
                            bands, Fnu = self.modelsdict[self.matched_parkeys_grid] 
                            matched_parkeys = [matched_parkeys]
                        else:
                            bands, Fnu = self.modelsdict[tuple(self.matched_parkeys_grid)] 
                        rest_bands = bands + np.log10((1+self.z))                         #Rest frame frequency
                        bandsf, Fnuf = f(rest_bands, Fnu, matched_parkeys[-1])  
                        bandsf = bandsf - np.log10((1+self.z))                            #Come back to frequency corrected by redshift

                        return bandsf, Fnuf

                    elif self.par_types[-2: ] == ['free', 'free'] and self.par_names[-2: ] == ['EBVbbb', 'alphaScat']: 
                        fcts=self.functions()
                        idxs= 0
                        f=fcts[self.functionidxs[idxs]]
                        bands, Fnu = self.modelsdict[tuple(self.matched_parkeys_grid)]
                        rest_bands = bands + np.log10((1+self.z))                         #Rest frame frequency
                        bandsf0, Fnuf0 = f(rest_bands[rest_bands < 16.685], Fnu[rest_bands < 16.685], matched_parkeys[-2])
                        bandsf =  np.concatenate((bandsf0, rest_bands[rest_bands >= 16.685])) - np.log10((1+self.z))   #Come back to redshifted frequency
                        Fnuf = np.concatenate((Fnuf0, Fnu[rest_bands >= 16.685]*10**(matched_parkeys[-1]/0.3838)))     #Add the effect of alpha_ox scatter

                        return bandsf, Fnuf

                    elif self.par_types[-3: ] == ['free', 'free', 'grid'] and self.par_names[-3: ] == ['EBVbbb', 'alphaScat', 'Gamma']: 
                        fcts=self.functions()
                        idxs= 0
                        f=fcts[self.functionidxs[idxs]]
                        bands, Fnu = self.modelsdict[tuple(self.matched_parkeys_grid)]
                        rest_bands = bands + np.log10((1+self.z))                         #Rest frame frequency
                        bandsf0, Fnuf0 = f(rest_bands[rest_bands < 16.685], Fnu[rest_bands < 16.685], matched_parkeys[-3])
                        bandsf =  np.concatenate((bandsf0, rest_bands[rest_bands >= 16.685])) - np.log10((1+self.z))   #Come back to redshifted frequency
                        Fnuf = np.concatenate((Fnuf0, Fnu[rest_bands >= 16.685]*10**(matched_parkeys[-2]/0.3838)))     #Add the effect of alpha_ox scatter

                        return bandsf, Fnuf
                    else: 
                        print('Error DICTIONARIES_AGNfitter.py: parameter type ',self.par_types, ' is unknown.')

    if MD.modelsettings['RADIO']== True:
        agnrad_parkeys = np.array(list(MD.AGN_RADFdict.keys()))
        galaxy_parnames, starburst_parnames,torus_parnames, bbb_parnames, agnrad_parnames, norm_parnames = MD.all_parnames
        galaxy_partypes, starburst_partypes,torus_partypes, bbb_partypes, agnrad_partypes, norm_partypes = MD.all_partypes 
        agnrad_obj=get_model(agnrad_parnames,agnrad_partypes,agnrad_parkeys,MD.AGN_RADFdict,MD.z, MD.AGN_RADfunctions ,model.AGN_RADfunctions)

    else:     
        galaxy_parnames, starburst_parnames,torus_parnames, bbb_parnames, norm_parnames = MD.all_parnames
        galaxy_partypes, starburst_partypes,torus_partypes, bbb_partypes, norm_partypes = MD.all_partypes 
        agnrad_obj = '-99.9'   #If there isn't AGN radio model create a false object, so the get model class always return 5 elements

    gal_obj =get_model(galaxy_parnames,galaxy_partypes,galaxy_parkeys, MD.GALAXYFdict, MD.z, MD.GALAXYfunctions, model.GALAXYfunctions)
    sb_obj =get_model(starburst_parnames,starburst_partypes,starburst_parkeys, MD.STARBURSTFdict,MD.z, MD.STARBURSTfunctions, model.STARBURSTfunctions)
    tor_obj=get_model(torus_parnames,torus_partypes,torus_parkeys, MD.TORUSFdict,MD.z, MD.TORUSfunctions, model.TORUSfunctions)
    bbb_obj=get_model(bbb_parnames,bbb_partypes,bbb_parkeys,MD.BBBFdict,MD.z, MD.BBBfunctions ,model.BBBfunctions)

    return gal_obj,sb_obj,tor_obj, bbb_obj, agnrad_obj


def dictkey_arrays_4plot(MD):
    """
    Summarizes the model dictionary keys and does the interpolation to nearest value in grid.
    used to be transporte to data

    ##input:

    ##output:
    """

    galaxy_parkeys= np.array(list(MD.GALAXYFdict.keys()))
    starburst_parkeys = np.array(list(MD.STARBURSTFdict.keys()))
    torus_parkeys = np.array(list(MD.TORUSFdict.keys()))
    bbb_parkeys = np.array(list(MD.BBBFdict.keys()))

#    class pick_obj:
    class get_model:
            def __init__(self, par_names, par_types, pars_modelkeys, modelsdict, z, functionidxs, functions):

                self.pars_modelkeys=pars_modelkeys.T
                self.pars_modelkeys_float =self.pars_modelkeys.astype(float)
                self.par_names = par_names
                self.par_types = par_types
                self.modelsdict = modelsdict
                self.functions = functions
                self.functionidxs=functionidxs
                self.z= z

            def pick_nD(self, pars_mcmc): 
                self.matched_parkeys = []
                self.matched_parkeys_grid = []

                if len(pars_mcmc)==1:
                    for i in range(len(pars_mcmc)):   
                        if self.par_types[i] == 'grid':
                            matched_idx =np.abs(self.pars_modelkeys_float-pars_mcmc[i]).argmin() #Choose the parameter value closest to that found by mcmc
                            matched_parkey = self.pars_modelkeys[matched_idx]
                            self.matched_parkeys =matched_parkey
                            self.matched_parkeys_grid=self.matched_parkeys
                        elif self.par_types[i] == 'free':
                            self.matched_parkeys=pars_mcmc[i]                                   #Values found by mcmc in case of free parameters
                            self.matched_parkeys_grid = self.pars_modelkeys[0]
                        else: 
                            print('Error DICTIONARIES_AGNfitter.py: parameter type ',self.par_types, ' is unknown.')

                else:
                    for i in range(len(pars_mcmc)):
                        if self.par_types[i] == 'grid': #if not 'grid'  
                            matched_idx =np.abs(self.pars_modelkeys_float[i]-pars_mcmc[i]).argmin() #Choose the parameter value closest to that found by mcmc
                            matched_parkey = self.pars_modelkeys[i][matched_idx]
                            self.matched_parkeys.append(matched_parkey)
                            self.matched_parkeys_grid.append(matched_parkey)
                        elif self.par_types[i] == 'free':
                            #print ('line 206 dic : Free partype')
                            self.matched_parkeys.append(pars_mcmc[i])                              #Values found by mcmc in case of free parameters
                            self.matched_parkeys_grid.append(self.pars_modelkeys[i, 0])
                        else:
                            print('Error DICTIONARIES_AGNfitter.py: parameter type ',self.par_types, ' is unknown.')

                    self.matched_parkeys=tuple(self.matched_parkeys)

            def get_fluxes(self,  matched_parkeys):    #From the dictionary of original models (without redshifting and filtering)
                    
                    if 'free' not in self.par_types:  
                        return self.modelsdict[matched_parkeys]


                    elif 'free' in self.par_types and self.par_names[-1] == 'EBVgal':   #This is for the case of EBVgal == free
                        fcts=self.functions()
                        idxs=0
                        f=fcts[self.functionidxs[idxs]]
                        rest_bands, Fnu = self.modelsdict[tuple(self.matched_parkeys_grid)] 
                        bandsf, Fnuf = f(10**rest_bands, Fnu, matched_parkeys[-1])    #Calzetti function need normal frequency (not log)

                        return np.log10(bandsf), Fnuf

                    elif 'free' in self.par_types and self.par_names[-1] == 'EBVbbb':  #This is for the case of EBVbbb == free
                        fcts=self.functions()
                        idxs=0
                        f=fcts[self.functionidxs[idxs]]
                        #R06 without X-Rays only have 1 parameter (EBV_bbb) and not a list of parameters so tuple() produce problems
                        if type(self.matched_parkeys_grid) != list:          
                            rest_bands, Fnu = self.modelsdict[self.matched_parkeys_grid] 
                            matched_parkeys = [matched_parkeys]
                        else:
                            rest_bands, Fnu = self.modelsdict[tuple(self.matched_parkeys_grid)] 
                        bandsf, Fnuf = f(rest_bands, Fnu, matched_parkeys[-1])    

                        return bandsf, Fnuf

                    elif self.par_types[-2: ] == ['free', 'free'] and self.par_names[-2: ] == ['EBVbbb', 'alphaScat']:
                        fcts=self.functions()
                        idxs=0
                        f=fcts[self.functionidxs[idxs]]
                        rest_bands, Fnu = self.modelsdict[tuple(self.matched_parkeys_grid)]
                        bandsf0, Fnuf0 = f(rest_bands[rest_bands < 16.685], Fnu[rest_bands < 16.685], matched_parkeys[-2])
                        bandsf =  np.concatenate((bandsf0, rest_bands[rest_bands >= 16.685])) 
                        Fnuf = np.concatenate((Fnuf0, Fnu[rest_bands >= 16.685]*10**(matched_parkeys[-1]/0.3838)))   #Add the effect of alpha_ox scatter
                        return bandsf, Fnuf

                    elif self.par_types[-3: ] == ['free', 'free', 'grid'] and self.par_names[-3: ] == ['EBVbbb', 'alphaScat', 'Gamma']: 
                        fcts=self.functions()
                        idxs=0
                        f=fcts[self.functionidxs[idxs]]
                        rest_bands, Fnu = self.modelsdict[tuple(self.matched_parkeys_grid)]
                        bandsf0, Fnuf0 = f(rest_bands[rest_bands < 16.685], Fnu[rest_bands < 16.685], matched_parkeys[-3])
                        bandsf =  np.concatenate((bandsf0, rest_bands[rest_bands >= 16.685])) 
                        Fnuf = np.concatenate((Fnuf0, Fnu[rest_bands >= 16.685]*10**(matched_parkeys[-2]/0.3838)))   #Add the effect of alpha_ox scatter
                        return bandsf, Fnuf

                    else: 
                        print('Error DICTIONARIES_AGNfitter.py: parameter type ',self.par_types, ' is unknown.')

    if MD.modelsettings['RADIO']== True:
        agnrad_parkeys = np.array(list(MD.AGN_RADFdict.keys()))
        galaxy_parnames, starburst_parnames,torus_parnames, bbb_parnames, agnrad_parnames, norm_parnames = MD.all_parnames
        galaxy_partypes, starburst_partypes,torus_partypes, bbb_partypes, agnrad_partypes, norm_partypes = MD.all_partypes 
        agnrad_obj=get_model(agnrad_parnames,agnrad_partypes,agnrad_parkeys,MD.AGN_RADFdict_4plot,MD.z, MD.AGN_RADfunctions, model.AGN_RADfunctions)

    else:     
        galaxy_parnames, starburst_parnames,torus_parnames, bbb_parnames, norm_parnames = MD.all_parnames
        galaxy_partypes, starburst_partypes,torus_partypes, bbb_partypes, norm_partypes = MD.all_partypes 
        agnrad_obj = '-99.9'   #If there isn't AGN radio model create a false object, so the get model class always return 5 elements

    gal_obj =get_model(galaxy_parnames,galaxy_partypes,galaxy_parkeys, MD.GALAXYFdict_4plot, MD.z, MD.GALAXYfunctions, model.GALAXYfunctions)
    sb_obj =get_model(starburst_parnames,starburst_partypes,starburst_parkeys, MD.STARBURSTFdict_4plot,MD.z, MD.STARBURSTfunctions, model.STARBURSTfunctions)
    tor_obj=get_model(torus_parnames,torus_partypes,torus_parkeys, MD.TORUSFdict_4plot,MD.z, MD.TORUSfunctions, model.TORUSfunctions)
    bbb_obj=get_model(bbb_parnames,bbb_partypes,bbb_parkeys,MD.BBBFdict_4plot,MD.z, MD.BBBfunctions ,model.BBBfunctions)

    return gal_obj,sb_obj,tor_obj, bbb_obj, agnrad_obj



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


