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
import matplotlib
matplotlib.use('Agg')
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
from functions.DICTIONARIES_AGNfitter import MODELSDICT
from astropy import units as u
from types import *

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


def MAKE_model_dictionary(cat, filters, models, clobbermodel=False):
    """
    Create model dictionary for all redshifts in z-array and filters
    input 
        cat - catalog settings
        filters - filter settings
        clobbermodel - remove any existing dictionary (default - False)
    ouput
        modelsdict
    """
    
    t0= time.time()
    modelsdict_name = cat['path']+models['path']+ filters['filterset']+models['modelset']

    if clobbermodel and os.path.lexists(modelsdict_name) : ## If overwriting of model dictionary mode is on
        print "> Overwriting model dictionary (-o) "+modelsdict_name
        print '_______________________'
        os.system('rm -rf '+ modelsdict_name)
  
    if not os.path.lexists(modelsdict_name): ## If model dictionary does not exist yet

        print '> Constructing new MODELS DICTIONARY:'
        print 'SAVED AS : ', modelsdict_name
        print 'FILTERSET USED : ', filters['filterset']
        print '________________________'
        print 'This process might take some time but you have to do it only once.'
        print 'If you interrupt it, please trash the empty file created.'
        print ''

        mydict = MODELSDICT( modelsdict_name, cat['path'], filters, models)
        mydict.build()
        f = open(mydict.filename, 'wb')
        cPickle.dump(mydict, f, protocol=2)
        f.close()

        print '________________________'
        print 'The models dictionary ' + modelsdict_name +' has been created.'\
              '(%.2g min elapsed)'% ((time.time() - t0)/60.)

    else: ## If model dictionary exists and you want to reuseit

        mydict = cPickle.load(file(modelsdict_name, 'rb'))

        print '________________________'
        print 'MODELS DICTIONARY currently in use:'
        print 'SAVED AS : ', mydict.filename
        print 'FILTERSET USED : ', mydict.filterset_name
        
        test_settingschanges= [mydict.filters_list[i]==filters[i]  for i in filters.keys() if  type(filters[i]) is BooleanType]

        if False in test_settingschanges : ##  compare nr of data bands with n of model filters
            print '________________________'
            print '*** WARNING ***'  
            print 'You have changed entries in your FILTER_settings without updating the filterset name.'
            print 'The MODEL DICTIONARY will NOT be updated correctly.'
            print ''
            py3 = sys.version_info[0] > 2 #creates boolean value for test that Python major version > 2
            if py3:
              response = input(">> You want to continue with the old model dictionary nevertheless? (yes/no)")
            else:
              response = raw_input(">> You want to continue with the old model dictionary nevertheless? (yes/no)")
            while True:
                if response== 'no':
                    sys.exit('\nTo update the MODEL DICTIONARY , you have these options: \n\
                        * change the filterset name (filters["filterset"]), or \n\
                        * run the code in model-overwriting mode (-o)')
                elif response== 'yes':
                    print '________________________'
                    print 'MODELS DICTIONARY currently in use:'
                    print 'SAVED AS : ', mydict.filename
                    print 'FILTERSET USED : ', mydict.filterset_name
                    break
                else:
                    print "Please answer (yes/no)"
                    break

    Modelsdict = mydict.MD ## MD is the model dictionary saved as a class of the method MODELSDICT

    bands_in_cat =len(cat['freq/wl_list'])
    bands_in_dict = len(Modelsdict[Modelsdict.keys()[0]][1][Modelsdict[Modelsdict.keys()[0]][1].keys()[0]][0])

    if bands_in_cat!= bands_in_dict :

        sys.exit('________________________\n*** ERROR ***\n'+ \
                'Number of bands in catalog ('+str(bands_in_cat)+'), does not match the numbers of filters in models (' + str(bands_in_dict)+')\n'+\
                'This will produce a mismatched fitting. Make sure the filterset contains only/all the photometric bands corresponding your catalog.\n'+ \
                'They do NOT need to be sorted in the same order.')
        
    return Modelsdict


def RUN_AGNfitter_onesource_independent( line, data_obj, filtersz, models, clobbermodel=False):
    """
    Main function for fitting a single source in line and create it's modelsdict independently.
    """
    
    mc = MCMC_settings()
    out = OUTPUT_settings()
    data = DATA(data_obj,line)

    print ''
    print '________________________'    
    print 'Fitting sources from catalog: ', data.catalog 
    print '- Sourceline: ', line
    print '- Sourcename: ', data.name


    ## 0. CONSTRUCT DICTIONARY for this redshift
    t0= time.time()

    ## needs a list/array of z
    filtersz['dict_zarray'] = [data.z]

    ## save the dictionary for this source in the OUTPUT folder for this source
    ## create this source output folder if it doesn't exist
    if not os.path.lexists(cat['output_folder']+str(data.name)):
        os.system('mkdir -p ' + cat['output_folder'] +str(data.name))
    dictz = cat['output_folder'] +str(data.name) +'/MODELSDICT_' + str(data.name) 
    ## remove this source modelsdict if it already exists and we want to remove it
    if clobbermodel and os.path.lexists(dictz):
        os.system('rm -rf '+dictz)
        print "removing source model dictionary "+dictz
    try:  
        if not os.path.lexists(dictz):
            zdict = MODELSDICT(dictz, cat['path'], filtersz, models)
            zdict.build()
            f = open(zdict.filename, 'wb')
            cPickle.dump(zdict, f, protocol=2)
            f.close()
            print '_____________________________________________________'
            print 'For this dictionary creation %.2g min elapsed'% ((time.time() - t0)/60.)
        else:
            zdict = cPickle.load(file(dictz, 'rb'))
        Modelsdictz = zdict.MD


        data.DICTS(filtersz, Modelsdictz)



        P = parspace.Pdict (data)   # Dictionary with all parameter space specifications.
                                    # From PARAMETERSPACE_AGNfitter.py

        t1= time.time()
        MCMC_AGNfitter.main(data, P, mc)    
        PLOTandWRITE_AGNfitter.main(data,  P,  out, models)

        # try:
        #     PLOTandWRITE_AGNfitter.main(data,  P,  out, models)
        # except:
        #     MCMC_AGNfitter.main(data, P, mc)    
        #     try:    
        #         PLOTandWRITE_AGNfitter.main(data,  P,  out, models)
        #     except:
        #         print line, 'is problematic to plot and write'

        print '_____________________________________________________'
        print 'For this fit %.2g min elapsed'% ((time.time() - t1)/60.)

    except EOFError: 
        
        print 'Line ',line,' cannot be fitted.'

def RUN_AGNfitter_onesource( line, data_obj, models):
    """
    Main function for fitting a single source in line 'line'.
    """
    
    mc = MCMC_settings()
    out = OUTPUT_settings()
    data = DATA(data_obj,line)
    try:
        data.DICTS(filters, Modelsdict)
                #PLOTandWRITE_AGNfitter.main(data,  P,  out, models)
    except EOFError:
        print 'Line ',line,' cannot be fitted.'

    P = parspace.Pdict (data)  # Dictionary with all parameter space especifications.
                                # From PARAMETERSPACE_AGNfitter.py

    print ''
    print '________________________'      
    print 'Fitting sources from catalog: ', data.catalog 
    print '- Sourceline: ', line
    print '- Sourcename: ', data.name

    t1= time.time()
    MCMC_AGNfitter.main(data, P, mc)        
    PLOTandWRITE_AGNfitter.main(data,  P,  out, models)

    # try:
    #     MCMC_AGNfitter.main(data, P, mc)        
    #     PLOTandWRITE_AGNfitter.main(data,  P,  out, models)
    #     print 'Done already'        
    # except:
    #     print 'Not done yet'
    #     MCMC_AGNfitter.main(data, P, mc)        
    #     PLOTandWRITE_AGNfitter.main(data,  P,  out, models)

    print '_____________________________________________________'
    print 'For this fit %.2g min elapsed'% ((time.time() - t1)/60.)
    return

    
def multi_run_wrapper(args):
    """
    wrapper to allow calling RUN_AGNfitter_onesource in pool.map
    """
    return RUN_AGNfitter_onesource(*args)
def multi_run_wrapper_indep(args):
    """
    wrapper to allow calling RUN_AGNfitter_onesource in pool.map
    """
    return RUN_AGNfitter_onesource_independent(*args)

def RUN_AGNfitter_multiprocessing(processors, data_obj, modelsdict, indep_bool=False, filters='None'):
    """
    Main function for fitting all sources in a large catalog.
    Splits the job of running the large number of sources
    into a chosen number of processors.
    """
    
    if indep_bool==False:
        nsources = data_obj.cat['nsources']
        
        print "processing all {0:d} sources with {1:d} cpus".format(nsources, processors)
        
        pool = mp.Pool(processes = processors)
        catalog_fitting = pool.map(multi_run_wrapper, itertools.izip(range(nsources), itertools.repeat(data_obj), itertools.repeat(modelsdict)))
        pool.close()
        pool.join()
    else:
        nsources = data_obj.cat['nsources']
        
        print "processing all {0:d} sources with {1:d} cpus".format(nsources, processors)
        
        pool = mp.Pool(processes = processors)
        catalog_fitting = pool.map(multi_run_wrapper_indep, itertools.izip(range(nsources), itertools.repeat(data_obj), itertools.repeat(filters),itertools.repeat(modelsdict)))
        pool.close()
        pool.join()

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
       filters= FILTERS_settings()
       models= MODELS_settings()
    except NameError:
        print "Something is wrong with your setting file"
        sys.exit(1)
        

    data_ALL = DATA_all(cat, filters)
    data_ALL.PROPS()

    ## make sure the output paths exist
    if not os.path.isdir(cat['output_folder']):
        os.system('mkdir -p '+os.path.abspath(cat['output_folder']))
    # abspath is needed because 'dict_path' is a file
    modelsdict_name = cat['path']+models['path']+ filters['filterset']+models['modelset']
    mpath = modelsdict_name.replace(os.path.basename(modelsdict_name),'')

    if not os.path.isdir(mpath):
        os.system('mkdir -p '+os.path.abspath(mpath))
    

    # run for one source only and construct dictionary only for this source
    if args.independent:
        if args.ncpu>1.:
            
            RUN_AGNfitter_multiprocessing(args.ncpu, data_ALL, models, indep_bool=True, filters=filters)

        elif args.sourcenumber >= 0:
            RUN_AGNfitter_onesource_independent(args.sourcenumber, data_ALL, filters, models, clobbermodel=clobbermodel)
        else:
            for i in range(0, 110, 1):
                RUN_AGNfitter_onesource_independent(i, data_ALL, filters, models, clobbermodel=clobbermodel)
            
        
    else:
        # make/read the model dictionary
        Modelsdict = MAKE_model_dictionary(cat, filters, models, clobbermodel=clobbermodel)

        # a single source is specified
        if args.sourcenumber >= 0 and args.independent==False:
            RUN_AGNfitter_onesource(args.sourcenumber, data_ALL, models)
        else:
            RUN_AGNfitter_multiprocessing(args.ncpu, data_ALL, models)
        
        
    print '======= : ======='
    print 'Process finished.'
