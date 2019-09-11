'''
AGNfitter setting file:

required:
  CATALOG_settings
  FILTERS_settings
  MCMC_settings
  OUTPUT_settings

For default use (test example with 2 redshifts and default filter set)

Change only the functions which state 
***USER INPUT NEEDED***.
'''

def CATALOG_settings():

    """==================================
    ***USER INPUT NEEDED***

    Set the right values to be able to read your catalog's format.
    FITS option is not available yet.
    =================================="""


    cat = dict()


    ##GENERAL
    cat['path'] ='/Users/Gabriela/Desktop/AGNfitter/'  #path to the AGNfitter code
    cat['filename'] = cat['path']+'data/mockGalSED.dat'
    cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. 
    cat['name'] = 0                 ## If ASCII: Column index (int) of source IDs
                                    ## If FITS : Column name (str). E.g. 'ID'
    cat['redshift'] = 2             ## If ASCII:  Column index(int) of redshift 
                                     ## If FITS : Column name (str). E.g. z'

    ##FREQUENCIES/WAVELENGTHS 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['freq/wl_list'] = np.arange(3, 16, 1).tolist() #64 wise                                 
                                        ## If ASCII: List of column indexes (int), 
                                        ##           corresponding to freq/wl.                                  
    #cat['freq/wl_suffix'] = '_wl'      ## If FITS: common ending to wavelength column names
    cat['use_central_wavelength'] = True # Option to use central wavelength if no wavelengths in table
    # Otherwise:
    cat['freq/wl_format'] = 'frequency' ## Gives your catalog *observed*
                                         ## 'frequency' or 'wavelength'?
    cat['freq/wl_unit'] = u.Hz    ## Astropy unit of freq or wavelength

    ##FLUXES 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['flux_unit'] = u.Jy            ## Astropy unit of *flux* (astropy-units)
    cat['flux_list'] = np.arange(16, 29, 1).tolist()        
                                        ## If ASCII: List of column indexes (int)
    #cat['flux_suffix'] = '_f'          ## If FITS: Common ending of all flux column names (str)    
    cat['fluxerr_list'] = np.arange(29, 42, 1).tolist() 
                                        ## If ASCII: List of column indexes (int)
    #cat['fluxerr_suffix'] = '_e'       ## If FITS: common ending to fluxerr column names (str)

    ##NON-DETECTIONS                                        
    cat['ndflag_bool'] = False          ## Does you catalog has columns with flags 1(0) for 
                                        ## detections (nondetections)? 
    cat['ndflag_list'] = 'list'         ## If ASCII: List of column indexes (int)
                                        ## If FITS: List of column names (str)    

    ## COSTUMIZED WORKING PATHS
    cat['workingpath'] = cat['path']  # Allows for a working path other than the AGNfitter code path.
                                      # Will include:
                                            # dictionary of models 
                                            # SETTINGS_AGNFitter.py file  
                                            # OUTPUT
                                      # Specially needed in order not to alter git original repository
                                      # and when using an external processor.
                                      # Default: cat['path'] (same as AGNfitter code path) 

    cat['use_central_wavelength'] = True # Option to use central wavelength if no wavelengths in table

    cat['output_folder'] = '/Users/Gabriela/Desktop/AGNfitter_OUTPUT/MOCK/' #if no special OUTPUT folder, leave default


    return cat


def FILTERS_settings():

    """==================================
    Set the photometric bands included in your catalog,
    in order to integrate the models over their response curves.
    =================================="""

    filters = dict()

    filters['dict_zarray'] = [8] The grid of redshifts needed to fit your catalog
    filters['path'] = 'models/FILTERS/' 
    filters['filterset'] = 'filterv2' ## 'filterset_default' (for the test case),
    

   #  filters['ALMA_band7'] = [True,26]    
   #  filters['SPIRE500'] = [True, 25]
   #  filters['SPIRE350'] = [True, 24]
   #  filters['SPIRE250'] = [True, 23]
   #  filters['PACS160'] = [True, 22]
   #  filters['PACS100'] = [True, 21]
   #  filters['MIPS70'] = [True, 20]
   #  filters['MIPS24'] = [True, 19]
   # # filters['WISE4'] = [True, 22]
   # # filters['WISE3'] = [True, 21]
   #  filters['IRAC4'] = [True, 18]
   #  filters['IRAC3'] = [True, 17]
   # # filters['WISE2'] = [True, 18]   
    filters['IRAC2'] = [True, 12]
    filters['IRAC1'] = [True, 11]
 #   filters['WISE1'] = [True, 15]
    # filters['wircam_Ks'] = [True, 14] 
    # filters['HAWKI_K'] = [True, 13]
    # filters['MUSYC_K'] = [True, 12]
    # filters['MUSYC_H'] = [True,11]
    filters['Ks_2mass'] = [True,10]
    filters['H_2mass'] = [True,9]
    filters['F160W'] = [True,8]
    filters['F140W'] = [True, 7]
    filters['F105W'] = [True,6]
    filters['z_SDSS'] = [True, 5]
    filters['i_SDSS'] = [True, 4]
    filters['r_SDSS'] = [True, 3]
    filters['g_SDSS'] = [True, 2]
    filters['u_SDSS'] = [True, 1]
    filters['GALEX_1500'] = [True, 0]                                              ## for the user's case: customize, eg. filtersv1
 
    filters['add_filters'] = False # If 'True' please add them below in ADD FILTERS

    """==================================
    ADD FILTERS (optional)
    =================================="""

    ADDfilters = dict()
    ADDfilters['names'] = ['MUSYC_J', 'HAWKI_J', 'TENIS_J', 'MUSYC_H','MUSYC_K', 'ALMA_band7']    ## (string/list of strings)User especified filter names. 
                                ## If name has been previously used, an error message will appear. 
    ADDfilters['filenames'] = ['models/FILTERS/MUSYC/ecdfs.J.filt.dat.txt', 'models/FILTERS/Paranal/Paranal_HAWKI.J.dat.txt', \
                                'models/FILTERS/CHFT/CFHT_Wircam.J.dat.txt', 'models/FILTERS/MUSYC/ecdfs.H.filt.dat.txt',\
                                'models/FILTERS/MUSYC/ecdfs.K.filt.dat.txt', 'models/FILTERS/ALMA/Alma7.txt'] ## (string/list of strings) File names of the new filters. 
                                ## File format: 2 columns of 1) freq/wavelength 2) Throughput. 
                                ## Path assumed is the cat['path'] especified above. 
                                ## Example: 'models/FILTERS/my_new_filter.txt'
    ADDfilters['freq/wl_format'] = ['wavelength'] * len(ADDfilters['names']) ## Info about the column 1 of your filter file.
                                                                             ## Options: 'wavelength' or 'frequency'.    
    ADDfilters['freq/wl_unit'] =  [u.Angstrom]* len(ADDfilters['names']) ## (Astropy Unit) Info about the column 1 
                                                                         ## of your filter file. 
    ADDfilters['description'] = ['description_dummy']* len(ADDfilters['names']) ## (Str) Any description the user wants to give 
                                                                                ##  to the filter to add.

    filters['add_filters_dict']= ADDfilters

    return filters

def MODELS_settings():

    """==================================
    Work in progress
    =================================="""

    models = dict()
    models['path'] = 'models/' 
    models['modelset'] = 'modelsv1'


    models['GALAXY'] = 'BC03'   ### Current options:
                                ### 'BC03' (Bruzual & Charlot 2003)

    models['STARBURST'] = 'S17' ### Current options:
                                ### 'DH02_CE01' (Dale & Helou 2002 + Chary & Elbaz 2001)
                                ### 'S07' (Schreiber et al. 2017 (submitted))
                                ### 'S17_newmodel' (Schreiber et al. 2017, less dust)

    models['BBB'] ='R06' ### Current options:
                         ### 'R06' (Richards et al. 2006) ## Needs 2 manual changes in PARAMETERSPACE_AGNfitter.py
                         ### 'SN12' (Slone&Netzer 2012)
                         ### 'D12_S' (Done et al. 2012) for Schwarzschild BH, with x-ray predictions
                         ### 'D12_K' (Done et al. 2012) for Kerr BH, with x-ray predictions

    models['TORUS'] ='S04' ### Current options:
                           ### 'S04' (Silva et al. 2004)
                           ###NK0


    return models

def MCMC_settings():

    """==================================
    Set your preferences for the MCMC sampling.
    =================================="""

    mc = dict()

    mc['Nwalkers'] = 100  ## number of walkers 
    mc['Nburnsets']= 1   ## number of burn-in sets
    mc['Nburn'] = 5000 ## length of each burn-in sets
    mc['Nmcmc'] = 10000  ## length of each burn-in sets
    mc['iprint'] = 1000 ## show progress in terminal in steps of this many samples

    return mc

def OUTPUT_settings():

    """==================================
    Set your preferences for the production of OUTPUT files. 
    =================================="""

    out = dict()

    out['plot_format'] = 'pdf'

    #CHAIN TRACES
    out['plot_tracesburn-in'] = False    
    out['plot_tracesmcmc'] = True

    #BASIC OUTPUT
    out['Nsample'] = 3000 ## out['Nsample'] * out['Nthinning'] <= out['Nmcmc']
    out['Nthinning'] = 10 ## This describes thinning of the chain to sample
    out['writepar_meanwitherrors'] = True ##Write output values for all parameters in a file.
    out['plot_posteriortriangle'] = True ##Plot triangle with all parameters' PDFs?
    out['chi_squared_value'] = True ##Write chi-squared value of best-fit SED, excludes upper limits

    #INTEGRATED LUMINOSITIES
    out['calc_intlum'] = True  
    out['realizations2int'] = 100 #This process is very time consuming.
                                #Around 100-1000 is recomendend for computational reasons.
                                #If you want to plot posterior triangles of 
                                #the integrated luminosities, should be > 1000.
    out['plot_posteriortrianglewithluminosities'] = False  # requires out['calc_intlum']=True
    out['save_posterior_luminosities'] = False

    #INTEGRATION RANGES
    out['intlum_models'] = ['sb','bbb', 'bbbdered', 'gal', 'tor','sb']  #leave 'sb' always 
                                                                        #as first element
    out['intlum_freqranges_unit'] = u.micron   #Astropy unit 
    out['intlum_freqranges'] = np.array([[8.,1000.],[0.1,1.],[0.1,1.],[5.,15.],[5.,15.],[5.,15.]])
    out['intlum_names'] = ['LIR(8-1000)','Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(5-15)', 'Ltor(5-15)','Lsb(5-15)']

    #SED PLOTTING
    out['realizations2plot'] = 10

    out['plotSEDrealizations'] = True

    return out
