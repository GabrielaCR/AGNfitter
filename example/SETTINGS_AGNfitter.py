'''
AGNfitter setting file:
required:
  CATALOG_settings
  FILTERS_settings
  MCMC_settings
  OUTPUT_settings

For default use (test example with 2 redshifts and default filter set)
change only the functions which state 
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
    cat['path'] ='/Users/USER/Codes/AGNfitter/'  #path to the AGNfitter code


    cat['filename'] = cat['path']+'data/catalog_example.txt'
    cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. 
                              ## FITS option not available yet. 

    cat['name'] = 0#'ID'            ## If ASCII: Column index (int) of source IDs
                                        ## If FITS: not yet
    cat['redshift'] = 1#'z'              ## If ASCII:  Column index(int) of redshift
                                        ## If FITS: not yet
 
    ##FREQUENCIES/WAVELENGTHS 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['freq/wl_list'] = np.arange(2,56,3).tolist()                                  
                                        ## If ASCII: List of column indexes (int), 
                                        ##           corresponding to freq/wl.                                  
    #cat['freq/wl_suffix'] = '_wl'      ## If FITS: common ending to wavelength column names
    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*
                                         ## 'frequency' or 'wavelength'?
    cat['freq/wl_unit'] = u.Angstrom       ## Astropy unit of freq or wavelength

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

    ##OTHERS (user input ***not*** needed)
    cat['outpath'] = '/Users/USER/my_AGNfitter/'  # allow for an output path outside of the AGNfitter code path
    cat['output_folder'] =  cat['outpath'] +'OUTPUT/'#if no special OUTPUT folder, leave default
    cat['dict_path'] = cat['outpath'] + 'models/MODELSDICT_default' 


    return cat


def FILTERS_settings():

    """==================================
    Set the photometric bands included in your catalog,
    in order to integrate the models over their response curves.
    =================================="""

    filters = dict()

    filters['dict_zarray'] =np.array([0.283, 1.58])  # The grid of redshifts needed to fit your catalog
    filters['Bandset'] = 'BANDSET_default' # OPTIONS: 
                                           # 'BANDSET_default' (for testing)
                                           # 'BANDSET_settings' (choosing relevant filters below, as given by your catalog)
                                           # if your filter is not included, go to DICTIONARIES_AGNfitter to add.

    filters['SPIRE500']= True
    filters['SPIRE350']= True
    filters['SPIRE250']= True
    filters['PACS160']=False
    filters['PACS100']=False

    filters['MIPS160']=False      
    filters['MIPS70']=False    
    filters['MIPS24']=True

    filters['IRAC4']=True       
    filters['IRAC3']=True
    filters['IRAC2']=True
    filters['IRAC1']=True

    filters['WISE4']=False
    filters['WISE3']=False
    filters['WISE2']=False
    filters['WISE1']=False

    filters['Ks_2mass']=True
    filters['H_2mass']=True
    filters['J_2mass']=True

    filters['H_VISTA']=False
    filters['J_VISTA']=False
    filters['K_VISTA']=False
    filters['Y_VISTA']=True

    filters['u_SDSS']=False  
    filters['g_SDSS']=False
    filters['r_SDSS']=False
    filters['i_SDSS']=False  
    filters['z_SDSS']=False

    filters['g_SUBARU']=False
    filters['r_SUBARU']=True
    filters['i_SUBARU']=True  
    filters['z_SUBARU']=True
    filters['B_SUBARU']=True
    filters['V_SUBARU']=False

    filters['u_CHFT']=True  
    filters['g_CHFT']=False
    filters['r_CHFT']=False
    filters['i_CHFT']=False  
    filters['z_CHFT']=False

    filters['GALEX_2500']=True
    filters['GALEX_1500']=False

    return filters

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
