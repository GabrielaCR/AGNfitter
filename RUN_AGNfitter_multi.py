
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
    cat['freq/wl_list'] = np.arange(2,56,3).tolist()                                  
                                        ## If ASCII: List of column indexes (int), 
                                        ##           corresponding to freq/wl.
                                        ## If FITS: not yet
    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*
                                         ## 'frequency' or 'wavelength'?
    cat['freq/wl_unit'] = u.Angstrom       ## Astropy unit of freq or wavelength

    ##FLUXES 
    cat['fluxerr_list'] = np.arange(4,58,3).tolist() 
                                        ## If ASCII: List of column indexes (int)
                                        ## If FITS: List of column names (str)
    cat['flux_unit'] = u.Jy             ## Astropy unit of *flux* (astropy-units)
    cat['flux_list'] = np.arange(3,57,3).tolist()        
                                        ## If ASCII: List of column indexes (int)
                                        ## If FITS: Common ending of all fluxes-columns' names .

    ##NON-DETECTIONS                                        
    cat['ndflag_bool'] = False          ## Does you catalog has columns with flags 1(0) for 
                                        ## detections (nondetections)? 
    cat['ndflag_list'] = 'list'         ## If ASCII: List of column indexes (int)
                                        ## If FITS: List of column names (str)    

    ##OTHERS (user input ***not*** needed)
    cat['output_folder'] =  cat['path'] +'OUTPUT/'#if no special OUTPUT folder, leave default
    cat['dict_path'] = cat['path'] + 'models/MODELSDICT_default' 


    return cat


def MCMC_settings():

    """==================================
    Set your preferences for the MCMC sampling.
    =================================="""

    mc = dict()

    mc['Bandset'] = 'BANDSET_default' 
    mc['Nwalkers'] = 100  ## number of walkers 
    mc['Nburnsets']= 2   ## number of burn-in sets
    mc['Nburn'] = 3000 ## length of each burn-in sets
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
    out['plot_tracesburn-in'] = True    
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
data_ALL = DATA_all(cat)
data_ALL.PROPS()


## 0. CONSTRUCT DICTIONARY (not needed if default is used)

if not os.path.lexists(cat['dict_path']):

    MODELFILES.construct(cat['path'])

    mydict = MODELSDICT(cat['dict_path'], cat['path'])
    mydict.build()

Modelsdict = cPickle.load(file(cat['dict_path'], 'rb')) 


def RUN_AGNfitter_onesource( line, data_obj=data_ALL, modelsdict= Modelsdict):
        """
        Main function for fitting a single source in line 'line'.
        """
        
        mc = MCMC_settings()
        out = OUTPUT_settings()

        data = DATA(data_obj,line)
        data.DICTS(mc, Modelsdict)

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


RUN()
