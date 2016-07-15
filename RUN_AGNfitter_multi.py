
#!/usr/bin/env python2.7

"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
AGNfitter
    
==============================
    
Fitting SEDs of AGN in a MCMC Approach
G.Calistro Rivera, E.Lusso, J.Hennawi, D. Hogg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
This is the main script.

For default use 
change only the functions which state 
***USER INPUT NEEDED***.

"""

#PYTHON IMPORTS
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import time
import shelve
import multiprocessing as mp 
import cPickle

#AGNfitter IMPORTS
from functions import  MCMC_AGNfitter, PLOTandWRITE_AGNfitter
import functions.PARAMETERSPACE_AGNfitter as parspace
from functions.DATA_AGNfitter import DATA, DATA_all
from functions.PLOTandWRITE_AGNfitter import CHAIN, FLUXES_ARRAYS
from functions.DICTIONARIES_AGNfitter import MODELSDICT
import functions.CONSTRUCT_modelobjects as MODELFILES
from astropy import units as u



def CATALOG_settings():

    """==================================
    ***USER INPUT NEEDED***

    Set the right values to be able to read your catalog's format.
    FITS option is not available yet.
    =================================="""


    cat = dict()


    ##GENERAL
    cat['path'] ='/Users/USER/AGNfitter/'  #YOUR OWN PATH


    cat['filename'] = 'data/catalog_example.txt'
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
    cat['output_folder'] =  cat['path'] +'OUTPUT/'#if no special OUTPUT folder, leave default
    cat['dict_path'] = cat['path'] + 'models/MODELSDICT_default' 


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


def RUN_AGNfitter_onesource( line, data_obj=data_ALL, modelsdict= Modelsdict):
        """
        Main function for fitting a single source in line 'line'.
        """
        
        mc = MCMC_settings()
        out = OUTPUT_settings()

        data = DATA(data_obj,line)
        data.DICTS(filters, Modelsdict)

        P = parspace.Pdict (data)  # Dictionary with all parameter space especifications.
                                   # From PARAMETERSPACE_AGNfitter.py

        print ''
        print 'Fitting sources from catalog: ', data.catalog 
        print '- Sourceline: ', line
        print '- Sourcename: ', data.name


        t1= time.time()

        MCMC_AGNfitter.main(data, P, mc)        
        PLOTandWRITE_AGNfitter.main(data,  P,  out)


        print '_____________________________________________________'
        print 'For this fit %.2g min elapsed'% ((time.time() - t1)/60.)


def RUN_AGNfitter_multiprocessing(processors,data_obj=data_ALL):
    """
    Main function for fitting all sources in a large catalog.
    Splits the job of running the large number of sources
    into a chosen number of processors.
    """
    cat = CATALOG_settings()
     
    with open(cat['filename'], 'r') as f:
        lines = f.readlines()
        catalog_lines = len([l for l in lines if l.strip(' \n') != '']) 

    pool = mp.Pool(processes = processors)
    catalog_fitting = pool.map(RUN_AGNfitter_onesource,range(catalog_lines))
    pool.close()
    pool.join()
    ##WRITE ALL RESULST IN ONE TABLE




def RUN():

    """==========================================

    ***USER INPUT NEEDED***

    CHOOSE between one of the two versions:

    1. fit single or few sources (input: [int] sourceline)
    or 
    2. fit total catalog (input: [int] nr. of processors )

    (Comment the option you are not using.)
    ==========================================="""

    RUN_AGNfitter_onesource(0)
    #RUN_AGNfitter_multiprocessing(1)

    print '======= : ======='
    print 'Process finished.'


def header():
    print '              '
    print '             XXXX'
    print '___________ XXX _________________________________________________'
    print '            XX      '
    print '            X     '
    print '            X                       AGNfitter                     '
    print '         __ X __                    ---------                ' 
    print '     /**\   |   /**\                                          '
    print '... (*** =  o  = ***) ...........................................'                                
    print '     \**/__ | __\**/                                     '
    print '            X              Fitting SEDs of AGN and Galaxies  '
    print '            X             in a MCMC Approach '
    print '           xx              (Calistro Rivera et al. 2016)    '   
    print '          xx               '            
    print '_______ xxx______________________________________________________'
    print '     xxxx'
    print ''




"""--------------------------------------------"""

cat = CATALOG_settings()
filters= FILTERS_settings()
data_ALL = DATA_all(cat)
data_ALL.PROPS()


## 0. CONSTRUCT DICTIONARY (not needed if default is used)

if not os.path.lexists(cat['dict_path']):

    MODELFILES.construct(cat['path'])

    mydict = MODELSDICT(cat['dict_path'], cat['path'], filters)
    mydict.build()

Modelsdict = cPickle.load(file(cat['dict_path'], 'rb')) 


if __name__ == "__main__":
    header()
    RUN()
