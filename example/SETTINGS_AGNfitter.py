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
    cat['path'] ='/Users/gcalistr/Documents/AGNfitter/'  #path to the AGNfitter code
    cat['filename'] = cat['path']+'data/catalog_example.txt'
    cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. 
    cat['name'] = 0                 ## If ASCII: Column index (int) of source IDs
                                    ## If FITS : Column name (str). E.g. 'ID'
    cat['redshift'] = 1             ## If ASCII:  Column index(int) of redshift 
                                     ## If FITS : Column name (str). E.g. z'

   ##FREQUENCIES/WAVELENGTHS 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['freq/wl_list'] = np.arange(5,48,2).tolist()                                  
                                        ## If ASCII: List of column indexes (int), 
                                        ##           corresponding to freq/wl.                                  
    #cat['freq/wl_suffix'] = '_wl'      ## If FITS: common ending to wavelength column names

    cat['use_central_wavelength'] = True # Option to use central wavelength if no wavelengths in table

    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*
                                         ## 'frequency' or 'wavelength'?
    cat['freq/wl_unit'] = u.Angstrom     ## Astropy unit of freq or wavelength

    ##FLUXES 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['flux_in_magAB'] = False # Option to calculate flux and flux_error from magnitude AB.
    cat['flux_unit'] = u.Jy * 1e-3            ## Astropy unit of *flux* (astropy-units)
    cat['flux_list'] = np.arange(2,61,2).tolist()        
                                        ## If ASCII: List of column indexes (int)
    #cat['flux_suffix'] = '_f'          ## If FITS: Common ending of all flux column names (str)    
    cat['fluxerr_list'] = np.arange(3,62,2).tolist() 
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
                                      
    cat['output_folder'] =  cat['workingpath'] +'OUTPUT/' #if no special OUTPUT folder, leave default



    return cat


def FILTERS_settings():

    """==================================
    Set the photometric bands included in your catalog,
    in order to integrate the models over their response curves.
    =================================="""

    filters = dict()

    filters['dict_zarray'] =np.array([0.283, 1.58])  # Deprecated. The grid of redshifts needed to fit your catalog
    filters['path'] = 'models/FILTERS/' 
    filters['filterset'] = 'example_30datapointa' ## 'filterset_default' (for the test case),
                                               ## for the user's case: customize, eg. filtersv1

  
    filters['SPIRE500'] = [True, 29]
    filters['SPIRE350'] = [True, 28]
    filters['SPIRE250'] = [True, 27]
    filters['PACS160'] = [True, 26]
    filters['PACS100'] = [True, 25]

    filters['WISE4'] = [True, 24]
    filters['WISE3'] = [True, 23]
    filters['WISE2'] = [True, 22]
    filters['WISE1'] = [True, 21]

    filters['K_UKIDSS'] = [True, 20]
    filters['H_UKIDSS'] = [True, 19]
    filters['J_UKIDSS'] = [True, 18]
    filters['Y_UKIDSS'] = [True, 17]

    filters['Ks_VISTA'] = [True, 16]
    filters['H_VISTA'] = [True, 15]
    filters['J_VISTA'] = [True, 14]
    filters['Y_VISTA'] = [True, 13]
    filters['Z_VISTA'] = [True, 12]

    filters['g_PS1'] = [True, 7]
    filters['r_PS1'] = [True, 8]
    filters['i_PS1'] = [True, 9]
    filters['z_PS1'] = [True, 10]
    filters['y_PS1_good'] = [True, 11]

    filters['z_SDSS'] = [True, 6]
    filters['i_SDSS'] = [True, 5]
    filters['r_SDSS'] = [True, 4]
    filters['g_SDSS'] = [True, 3]
    filters['u_SDSS'] = [True, 2]
    filters['GALEX_2500'] = [True, 1]
    filters['GALEX_1500'] = [True, 0]


    filters['add_filters']= False # If 'True' please add them below in ADD FILTERS

    """==================================
    ADD FILTERS (optional)
    =================================="""

    ADDfilters=dict()
    ADDfilters['names'] = []    ## (string/list of strings)User especified filter names. 
                                ## If name has been previously used, an error message will appear. 
    ADDfilters['filenames'] =[] ## (string/list of strings) File names of the new filters. 
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


    models['GALAXY'] = 'BC03_metal'   ### Current options:
                                ### 'BC03' (Bruzual & Charlot 2003)
                                ### 'BC03_metal' (Bruzual & Charlot 2003), with metallicities
    models['STARBURST'] = 'S17_newmodel' ### Current options:
                                ### 'DH02_CE01' (Dale & Helou 2002 + Chary & Elbaz 2001)
                                ### 'S07' (Schreiber et al. 2017 (submitted))

    models['BBB'] ='R06' ### Current options:
                         ### 'R06' (Richards et al. 2006) ## Needs 2 manual changes in PARAMETERSPACE_AGNfitter.py
                         ### 'SN12' (Slone&Netzer 2012)
                         ### 'D12_S' (Done et al. 2012) for Schwarzschild BH, with x-ray predictions
                         ### 'D12_K' (Done et al. 2012) for Kerr BH, with x-ray predictions

    models['TORUS'] ='S04' ### Current options:
                           ### 'S04' (Silva et al. 2004)

    models['XRAYS'] = False ### If X-ray data is available and informative for the fit

    models['RADIO'] = False ### If radio data is available and informative for the fit

    models['PRIOR_energy_balance'] = True ### Default:True
                                          ### True: Sets a lower limit to the dust emission luminosity ('starburst' model)
                                          ### as given by the observed attenuation in the stellar component SED.
    models['PRIOR_AGNfraction'] = True  ### Default: True
                                        ### True: - *IF* blue/UV bands (around 1500 Angstrom) are 10 times higher than expected by the galaxy luminosity function by Parsa, Dunlop et al. 2014. 
                                        ###         this option rejects AGN-to-GAL ratios lower than 1 (log =0). It then applies a Gaussian prior probability with log ratio=2, with a sigma of 2.
                                        ###       - In this cases it also applies a Gaussian prior on the galaxy normalization, i.e. stellar mass (usually unconstrained in these cases) to 
                                        ###         populate physically expected ranges for QSO hosts -> 10^9 - 10^11. 
                                        ###       - *ELSE IF* blue/UV bands (around 1500 Angstrom) are below 10 times the expected value by Parsa, Dunlop et al. 2014. 
                                        ###         this option gives preference to galaxy contribution in the optical UV, with Gaussian prior probability centered on AGN to GALAXY log ratios of -1. 
                                        ###          and sigma 1, i.e. accretion disk is disfavoured at least the data strongly prefers it.
                                        ### False:- Non-informative prior
    models['PRIOR_galaxy_only'] = False ### Default:False 
                                        ### True: sets all AGN contribution to 0.ß
    return models

def MCMC_settings():

    """==================================
    Set your preferences for the MCMC sampling.
    =================================="""

    mc = dict()

    mc['Nwalkers'] = 100  ## number of walkers 
    mc['Nburnsets']= 2   ## number of burn-in sets
    mc['Nburn'] = 4000 ## length of each burn-in sets
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
    out['Nsample'] = 1000 ## out['Nsample'] * out['Nthinning'] <= out['Nmcmc']
    out['Nthinning'] = 10 ## This describes thinning of the chain to sample
    out['writepar_meanwitherrors'] = True ##Write output values for all parameters in a file.
    out['plot_posteriortriangle'] = False ##Plot triangle with all parameters' PDFs?

    #INTEGRATED LUMINOSITIES
    out['calc_intlum'] = True  
    out['save_posterior_luminosities']= False
    out['realizations2int'] = 100 #This process is very time consuming.
                                #Around 100-1000 is recomendend for computational reasons.
                                #If you want to plot posterior triangles of 
                                #the integrated luminosities, should be > 1000.
    out['plot_posteriortrianglewithluminosities'] = False  # requires out['calc_intlum']=True

    #INTEGRATION RANGES
    out['intlum_models'] = ['sb','bbb', 'bbbdered', 'gal', 'tor','sb']  #leave 'sb' always 
                                                                        #as first element
    out['intlum_freqranges_unit'] = u.micron   #Astropy unit 
    out['intlum_freqranges'] = np.array([[8.,1000.],[0.1,1.],[0.1,1.],[0.1,1.],[1.,30.],[1.,30.]])
    out['intlum_names'] = ['LIR(8-1000)','Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)','Lsb(1-30)']

    #SED PLOTTING
    out['realizations2plot'] = 10

    out['plotSEDrealizations'] = True

    return out
