"""%%%%%%%%%%%%%%%%%

        DATA_AGNFitter.py

%%%%%%%%%%%%%%%%%%

This script contains the class DATA, 
which administrate the catalog properties given by the user props(). 
It also helps transporting
the main information on the dictionaries (DICTS).
"""
import sys,os
import numpy as np
from math import exp,log,pi, sqrt
import matplotlib.pyplot as plt
from numpy import random,argsort,sqrt
import time
from scipy.integrate import quad, trapz
from astropy import constants as const
from astropy import units as u
from astropy.table import Table
from astropy.io import fits
import cPickle
import functions.MODEL_AGNfitter as model
import functions.DICTIONARIES_AGNfitter as dicts



class DATA_all:

    """
    Class DATA_all
    ---------------
    Object with data info for the total catalog.
    It reads and processes all information about the catalog.
    It returns arrays with all important values (sourcenames, redshift, etc)
    and gives it to the class DATA, which administrates it for each sourceline.
    
    input: catalogname
    bugs: Not ready to read FITS yet.

    """

    def __init__(self, cat):
        self.cat = cat
        #self.sourceline = sourceline
        self.catalog = cat['filename']
        if not os.path.lexists(cat['filename']):
            print 'ERROR: Catalog does not exist under this name '+cat['filename']
            sys.exit(1)
        self.path = cat['path']
        self.dict_path = cat['dict_path']
        self.output_folder = cat['output_folder']

    def PROPS(self):

        if self.cat['filetype'] == 'ASCII': 
            #read all columns
            column = np.loadtxt(self.catalog, skiprows=1, unpack=True)

            #properties
            self.name = column[self.cat['name']].astype(int)
            self.z = column[self.cat['redshift']].astype(float)
            self.dlum = np.array([model.z2Dlum(z) for z in self.z])

            #read all wavelengths, fluxes, fluerrors, flags
            freq_wl_cat_ALL = \
                np.array([column[c] for c in self.cat['freq/wl_list']])* self.cat['freq/wl_unit'] 
            flux_cat_ALL =\
                np.array([ca for ca in  column[self.cat['flux_list']] ])*self.cat['flux_unit']
            fluxerr_cat_ALL = \
                np.array([ce for ce in column[self.cat['fluxerr_list']]])*self.cat['flux_unit']
            if self.cat['ndflag_bool'] == True: 
                ndflag_cat_ALL = np.array(column[self.cat['ndflag_list']])

            nus_l=[]
            fluxes_l=[]
            fluxerrs_l=[]
            ndflag_l=[]

            nrBANDS, nrSOURCES= np.shape(flux_cat_ALL)
            
            self.cat['nsources'] = nrSOURCES

            ##Convert to right units but give back just values
            for j in range(nrSOURCES):
            
                freq_wl_cat= freq_wl_cat_ALL[:,j]
                flux_cat= flux_cat_ALL[:,j]
                fluxerr_cat= fluxerr_cat_ALL[:,j]

                if self.cat['freq/wl_format']== 'frequency' :
                    nus0 = np.log10(freq_wl_cat.to(u.Hz).value)
                if self.cat['freq/wl_format']== 'wavelength' :
                    nus0 = np.log10(freq_wl_cat.to(u.Hz, equivalencies=u.spectral()).value)

                fluxes0 = np.array(flux_cat.to(u.erg/ u.s/ (u.cm)**2 / u.Hz).value)
                fluxerrs0 = np.array(fluxerr_cat.to(u.erg/ u.s/(u.cm)**2/u.Hz).value)

                ## If columns with flags exist
                if self.cat['ndflag_bool'] == True: 
                    ndflag_cat0 = ndflag_cat_ALL[:,j]     
                    # If fluxerrs0 are not given (-99), we assume flux is an upper limit for a non detection.
                    # Upper limit flux is then represented for the fitting
                    # with a data point at uppflux/2, and an error of +- uppflux/2
                    # implying an uncertanty that ranges from [0,uppflux]  
                    ndflag_cat0[fluxerr_cat.value<=-99]= 0.                                 
                    fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]= \
                                                            fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]*0.5
                    fluxerrs0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]= \
                                                            fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]

                ## If NO columns with flags exist
                elif self.cat['ndflag_bool'] == False:
                    ndflag_cat0 = np.ones(np.shape(fluxes0))
                    # If fluxerrs0 are not given (-99), we assume flux is an upper limit for a non detection.
                    # Upper limit flux is then represented for the fitting
                    # with a data point at uppflux/2, and an error of +- uppflux/2
                    # implying an uncertanty that ranges from [0,uppflux]
                    ndflag_cat0[fluxerr_cat.value<=-99]= 0.
                    fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]=\
                                                             fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]*0.5
                    fluxerrs0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]= \
                                                            fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]
                    # If neither fluxes and fluxerrs are given (both -99), 
                    # these are considered as a non existant data point.

                ## Sort in order of frequency
                nus_l.append(nus0[nus0.argsort()])
                fluxes_l.append(fluxes0[nus0.argsort()])
                fluxerrs_l.append(fluxerrs0[nus0.argsort()])
                ndflag_l.append(ndflag_cat0[nus0.argsort()])


            self.nus = np.array(nus_l)
            self.fluxes = np.array(fluxes_l)
            self.fluxerrs = np.array(fluxerrs_l)
            self.ndflag = np.array(ndflag_l)

        elif self.cat['filetype'] == 'FITS': 

            #read all columns
            fitstable = Table.read(self.catalog)

            #properties
            self.name = fitstable[self.cat['name']].astype(int)
            self.z = fitstable[self.cat['redshift']].astype(float)
            self.dlum = np.array([model.z2Dlum(z) for z in self.z])

            #read all wavelengths, fluxes, fluerrors, flags
            colnames = fitstable.dtype.names
            # handle the case when the columns are strangley ordered in the fits file (i.e. not band1_wl, band1_f, band1_e, band2_wl, band1_f, band1_f, band2_e, etc)
            # if only their suffixes are different, sorting them should put them in the same order
            wl_cols = [ c for c in colnames if self.cat['freq/wl_suffix'] in c]
            flux_cols = [ w.replace(self.cat['freq/wl_suffix'], self.cat['flux_suffix']) for w in wl_cols ]
            flux_err_cols = [ w.replace(self.cat['freq/wl_suffix'], self.cat['fluxerr_suffix']) for w in wl_cols ]
            
            # check that the flux and error columns exist in the fits table
            # stop running if they don't
            if np.any(np.array([f not in colnames for f in flux_cols])):
                print 'wavelength columns exist without corresponding flux columns in fits file:'
                for f in flux_cols:
                    if f not  in colnames: print f
                sys.exit(1)
            if np.any(np.array([f not in colnames for f in flux_err_cols])):
                print 'wavelength columns exist without corresponding flux err columns in fits file:'
                for f in flux_err_cols:
                    if f not  in colnames: print f
                sys.exit(1)

            freq_wl_cat_ALL = \
                np.array([fitstable[c] for c in wl_cols])* self.cat['freq/wl_unit'] 
            flux_cat_ALL =\
                np.array([fitstable[ca] for ca in  flux_cols ])*self.cat['flux_unit']
            fluxerr_cat_ALL = \
                np.array([fitstable[ce] for ce in flux_err_cols ])*self.cat['flux_unit']
            if self.cat['ndflag_bool'] == True: 
                ndflag_cat_ALL = np.array(fitstable[self.cat['ndflag_list']])

            nus_l=[]
            fluxes_l=[]
            fluxerrs_l=[]
            ndflag_l=[]

            nrBANDS, nrSOURCES= np.shape(flux_cat_ALL)

            self.cat['nsources'] = nrSOURCES
            
            ##Convert to right units but give back just values
            for j in range(nrSOURCES):
            
                freq_wl_cat= freq_wl_cat_ALL[:,j]
                flux_cat= flux_cat_ALL[:,j]
                fluxerr_cat= fluxerr_cat_ALL[:,j]

                if self.cat['freq/wl_format']== 'frequency' :
                    nus0 = np.log10(freq_wl_cat.to(u.Hz).value)
                if self.cat['freq/wl_format']== 'wavelength' :
                    nus0 = np.log10(freq_wl_cat.to(u.Hz, equivalencies=u.spectral()).value)

                fluxes0 = np.array(flux_cat.to(u.erg/ u.s/ (u.cm)**2 / u.Hz).value)
                fluxerrs0 = np.array(fluxerr_cat.to(u.erg/ u.s/(u.cm)**2/u.Hz).value)

                ## If columns with flags exist
                if self.cat['ndflag_bool'] == True: 
                    ndflag_cat0 = ndflag_cat_ALL[:,j]     
                    # If fluxerrs0 are not given (-99), we assume flux is an upper limit for a non detection.
                    # Upper limit flux is then represented for the fitting
                    # with a data point at uppflux/2, and an error of +- uppflux/2
                    # implying an uncertanty that ranges from [0,uppflux]  
                    ndflag_cat0[flux_cat.value<=-99]= 0.                                 
                    fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]=\
                                                     fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]*0.5
                    fluxerrs0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]=\
                                                     fluxes0 [(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]

                ## If NO columns with flags exist
                elif self.cat['ndflag_bool'] == False:

                    ndflag_cat0 = np.ones(np.shape(fluxes0))
                    # If fluxerrs0 are not given (-99), we assume flux is an upper limit for a non detection.
                    # Upper limit flux is then represented for the fitting
                    # with a data point at uppflux/2, and an error of +- uppflux/2
                    # implying an uncertanty that ranges from [0,uppflux]  
                    ndflag_cat0[flux_cat.value<=-99]= 0.
                    fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]=\
                                                     fluxes0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]*0.5
                    fluxerrs0[(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]=\
                                                     fluxes0 [(fluxerr_cat.value<=-99)&(flux_cat.value>-99)]
                    # If neither fluxes and fluxerrs are given (both -99), 
                    # these are considered as a non existant data point.

                ## Sort in order of frequency
                nus_l.append(nus0[nus0.argsort()])
                fluxes_l.append(fluxes0[nus0.argsort()])
                fluxerrs_l.append(fluxerrs0[nus0.argsort()])
                ndflag_l.append(ndflag_cat0[nus0.argsort()])


            self.nus = np.array(nus_l)
            self.fluxes = np.array(fluxes_l)
            self.fluxerrs = np.array(fluxerrs_l)
            self.ndflag = np.array(ndflag_l)

class DATA():

    """
    Class DATA
    ----------
    Object with data info for once source.
    It recieves the catalog information obtained from 
    object from class DATA_all and administrates it for each sourceline.

    input: object of class DATA_all, sourceline
    bugs: Not ready to read FITS yet.

    """

    def __init__(self, data_all, line):

        catalog = data_all
        self.nus = catalog.nus[line]
        self.fluxes = catalog.fluxes[line]
        self.fluxerrs = catalog.fluxerrs[line]
        self.ndflag = catalog.ndflag[line]
        self.name = catalog.name[line]
        self.z =catalog.z[line]
        self.dlum = catalog.dlum[line]

        self.cat = catalog.cat
        #self.sourceline = sourceline
        self.catalog = catalog.cat['filename']
        if not os.path.lexists(catalog.cat['filename']):
            print 'Catalog does not exist under this name.'
        self.path = catalog.cat['path']
        self.dict_path = catalog.cat['dict_path']
        self.output_folder = catalog.cat['output_folder']


    def DICTS(self, filters, Modelsdict):
        """
        Helps transporting the dictionary content
        corresponding to the redshift of the source
        """

        z_array = np.array(list(Modelsdict.keys()))
        idx = (np.abs(z_array.astype(float)-self.z)).argmin()
        z_key = z_array[idx] 

        self.filterdict = dicts.filter_dictionaries(filters['Bandset'], self.path, filters)   
        self.dict_modelfluxes = Modelsdict[z_key]
        self.dictkey_arrays = dicts.dictkey_arrays(self.dict_modelfluxes)
        
        print 'Filter set contains {:d} bands'.format(len(self.filterdict[0]))
        m = Modelsdict[z_key]
        bands = m[0][m[0].keys()[0]][0]
        print 'Model sets contains {:d} bands'.format(len(bands))
        
        
