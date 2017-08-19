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
    cat['path'] ='/Users/USER/AGNfitter/'  #path to the AGNfitter code
    cat['filename'] = cat['path']+'data/catalog_example.txt'
    cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. 
    cat['name'] = 0                 ## If ASCII: Column index (int) of source IDs
                                    ## If FITS : Column name (str). E.g. 'ID'
    cat['redshift'] = 1             ## If ASCII:  Column index(int) of redshift 
                                     ## If FITS : Column name (str). E.g. z'

    ##FREQUENCIES/WAVELENGTHS 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['freq/wl_list'] = np.arange(2,56,3).tolist()                                  
                                        ## If ASCII: List of column indexes (int), 
                                        ##           corresponding to freq/wl.                                  
    #cat['freq/wl_suffix'] = '_wl'      ## If FITS: common ending to wavelength column names
    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*
                                         ## 'frequency' or 'wavelength'?
    cat['freq/wl_unit'] = u.Angstrom     ## Astropy unit of freq or wavelength

    ##FLUXES 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['flux_unit'] = u.Jy             ## Astropy unit of *flux* (astropy-units)
    cat['flux_list'] = np.arange(3,57,3).tolist()        
                                        ## If ASCII: List of column indexes (int)
    #cat['flux_suffix'] = '_f'          ## If FITS: Common ending of all flux column names (str)    
    cat['fluxerr_list'] = np.arange(4,58,3).tolist() 
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

    filters['dict_zarray'] =np.array([0.283, 1.58])  # The grid of redshifts needed to fit your catalog
    filters['path'] = 'models/FILTERS/' 
    filters['filterset'] = 'filterset_default' ## 'filterset_default' (for the test case),
                                               ## for the user's case: customize, eg. filtersv1
    
    filters['SPIRE500'] = True
    filters['SPIRE350'] = True
    filters['SPIRE250'] = True
    filters['PACS160'] = True
    filters['PACS100'] = True
    filters['MIPS160'] = True
    filters['MIPS70'] = True
    filters['MIPS24'] = True
    filters['IRAC4'] = True
    filters['IRAC3'] = True
    filters['IRAC2'] = True
    filters['IRAC1'] = True
    filters['WISE4'] = True
    filters['WISE3'] = True
    filters['WISE2'] = True
    filters['WISE1'] = True
    filters['Ks_2mass'] = True
    filters['H_2mass'] = True
    filters['J_2mass'] = True
    filters['H_VISTA'] = True
    filters['J_VISTA'] = True
    filters['K_VISTA'] = True
    filters['Y_VISTA'] = True
    filters['u_CHFT'] = True
    filters['g_CHFT'] = True
    filters['r_CHFT'] = True
    filters['i_CHFT'] = True
    filters['z_CHFT'] = True
    filters['u_SDSS'] = True
    filters['g_SDSS'] = True
    filters['r_SDSS'] = True
    filters['i_SDSS'] = True
    filters['z_SDSS'] = True
    filters['g_SUBARU'] = True
    filters['r_SUBARU'] = True
    filters['i_SUBARU'] = True
    filters['z_SUBARU'] = True
    filters['B_SUBARU'] = True
    filters['V_SUBARU'] = True
    filters['GALEX_2500'] = True
    filters['GALEX_1500'] = True
    filters['MUSYC_U'] = True
    filters['MUSYC_B'] = True
    filters['MUSYC_V'] = True
    filters['MUSYC_R'] = True
    filters['MUSYC_I'] = True
    filters['MUSYC_z'] = True

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


    models['GALAXY'] = 'BC03'   ### Current options:
                                ### 'BC03' (Bruzual & Charlot 2003)

    models['STARBURST'] = 'S17' ### Current options:
                                ### 'DH02_CE01' (Dale & Helou 2002 + Chary & Elbaz 2001)
                                ### 'S07' (Schreiber et al. 2017 (submitted))

    models['BBB'] ='R06' ### Current options:
                         ### 'R06' (Richards et al. 2006)

    models['TORUS'] ='S04' ### Current options:
                           ### 'S04' (Silva et al. 2004)


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
