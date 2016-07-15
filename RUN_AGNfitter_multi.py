#!/usr/bin/env python2.7

"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
AGNfitter
    
==============================
    
Fitting SEDs of AGN in a MCMC Approach
G.Calistro Rivera, E.Lusso, J.Hennawi, D. Hogg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
This is the main script.


"""

#PYTHON IMPORTS
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import time
import shelve
import multiprocessing as mp 
import itertools
import cPickle
import argparse

#AGNfitter IMPORTS
from functions import  MCMC_AGNfitter, PLOTandWRITE_AGNfitter
import functions.PARAMETERSPACE_AGNfitter as parspace
from functions.DATA_AGNfitter import DATA, DATA_all
from functions.PLOTandWRITE_AGNfitter import CHAIN, FLUXES_ARRAYS
from functions.DICTIONARIES_AGNfitter import MODELSDICT
import functions.CONSTRUCT_modelobjects as MODELFILES
from astropy import units as u


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
    return


def MAKE_model_dictionary(cat, filters, clobbermodel=False):
    ## 0. CONSTRUCT DICTIONARY (not needed if default is used)
    t0= time.time()

    if clobbermodel and os.path.lexists(cat['dict_path']):
        print "removing model dictionary "+cat['dict_path']
        os.system('rm -rf '+ cat['dict_path'])
        
    if not os.path.lexists(cat['dict_path']):

        MODELFILES.construct(cat['path'])

        mydict = MODELSDICT(cat['dict_path'], cat['path'], filters)
        mydict.build()
        
        print '_____________________________________________________'
        print 'For this dictionary creation %.2g min elapsed'% ((time.time() - t0)/60.)

    Modelsdict = cPickle.load(file(cat['dict_path'], 'rb'))
    
    return Modelsdict


def RUN_AGNfitter_onesource_independent( line, data_obj, filtersz, clobbermodel=False):
    """
    Main function for fitting a single source in line 'line' and create it's modelsdict independently.
    """
    
    mc = MCMC_settings()
    out = OUTPUT_settings()

    data = DATA(data_obj,line)

    print ''
    print 'Fitting sources from catalog: ', data.catalog 
    print '- Sourceline: ', line
    print '- Sourcename: ', data.name


    ## 0. CONSTRUCT DICTIONARY for this redshift
    t0= time.time()

    # needs a list/array of z
    filtersz['dict_zarray'] = [data.z]

    # add a suffix for this source dictionary
    dictz = cat['dict_path'] + '_' + str(data.name) 
    if clobbermodel and os.path.lexists(dictz):
        os.system('rm -rf '+dictz)
        print "removing source model dictionary "+dictz
      
    if not os.path.lexists(dictz):
        zdict = MODELSDICT(dictz, cat['path'], filtersz)
        zdict.build()
        print '_____________________________________________________'
        print 'For this dictionary creation %.2g min elapsed'% ((time.time() - t0)/60.)

    Modelsdictz = cPickle.load(file(dictz, 'rb')) 

    data.DICTS(filtersz, Modelsdictz)


    P = parspace.Pdict (data)   # Dictionary with all parameter space especifications.
                                # From PARAMETERSPACE_AGNfitter.py

    t1= time.time()

    MCMC_AGNfitter.main(data, P, mc)        
    PLOTandWRITE_AGNfitter.main(data,  P,  out)


    print '_____________________________________________________'
    print 'For this fit %.2g min elapsed'% ((time.time() - t1)/60.)
    return

def RUN_AGNfitter_onesource( line, data_obj, modelsdict):
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
    return

    

def RUN_AGNfitter_multiprocessing(processors, data_obj, modelsdict):
    """
    Main function for fitting all sources in a large catalog.
    Splits the job of running the large number of sources
    into a chosen number of processors.
    """
    
    def multi_run_wrapper(args):
        return RUN_AGNfitter_onesource(*args)
    
    
    nsources = data_obj.cat['nsources']
    
    print "processing all {0:d} sources with {1:d} cpus".format(nsources, processors)
    
    pool = mp.Pool(processes = processors)
    catalog_fitting = pool.map(multi_run_wrapper, itertools.izip(range(nsources), itertools.repeat(data_obj), itertools.repeat(modelsdict)))
    pool.close()
    pool.join()
    ##WRITE ALL RESULST IN ONE TABLE
    return


if __name__ == "__main__":
    header()
  
    parser = argparse.ArgumentParser()
    
    parser.add_argument("AGNfitterSettings", type=str, help="AGNfitter settings file")
    parser.add_argument("-c","--ncpu", type=int, default=1, help="number of cpus to use for multiprocessing")
    parser.add_argument("-n", "--sourcenumber", type=int, default=-1, help="specify a single source number to run (this is the line number in hte catalogue not the source id/name)")
    parser.add_argument("-i","--independent", action="store_true", help="run independently per source, i.e. do not create a global model dictionary")
    parser.add_argument("-o","--overwrite", action="store_true", help="overwrite model files")
    
    
    
    args = parser.parse_args()
    
    execfile(args.AGNfitterSettings)
    
    if args.overwrite:
      clobbermodel = True
    else:
      clobbermodel = False
    
    try:
	cat = CATALOG_settings()
    except NameError:
        print "Something is wrong with your setting file"
        sys.exit(1)
        
    filters= FILTERS_settings()
    data_ALL = DATA_all(cat)
    data_ALL.PROPS()

    ## make sure the output paths exist
    if not os.path.isdir(cat['output_folder']):
        os.system('mkdir -p '+os.path.abspath(cat['output_folder']))
    # abspath is needed because 'dict_path' is a file
    mpath = cat['dict_path'].replace(os.path.basename(cat['dict_path']),'')
    if not os.path.isdir(mpath):
        os.system('mkdir -p '+os.path.abspath(mpath))
    

    # run for once source only and construct dictionary only for this source
    if args.independent:
        RUN_AGNfitter_onesource_independent(args.sourcenumber, data_ALL, filters, clobbermodel=clobbermodel)
        
        
    else:
        # make/read the model dictionary
        Modelsdict = MAKE_model_dictionary(cat, filters, clobbermodel=clobbermodel)

        # a single source is specified
        if args.sourcenumber >= 0:
            RUN_AGNfitter_onesource(args.sourcenumber, data_ALL, Modelsdict)
        else:
            RUN_AGNfitter_multiprocessing(args.ncpu, data_ALL, Modelsdict)
        
        
    print '======= : ======='
    print 'Process finished.'
