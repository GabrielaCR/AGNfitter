
"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         MOCKSEDs_AGNfitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains all functions which are needed to construct the total model of AGN. 
The functions here translate the parameter space points into total fluxes dependin on the models chosen.

Functions contained here are the following:

STARBURST_nf
BBB_nf
GALAXY_nf
TORUS_nf

"""

import numpy as np
from math import exp,pi, sqrt
import matplotlib.pyplot as plt
import time
import cPickle
from astropy.table import Table
from astropy.io import fits, ascii
from scipy.integrate  import quad, trapz
import astropy.constants as const
import astropy.units as u
import itertools

import MODEL_AGNfitter as model


		## After 83
        ##########################################erase!!!!!
        i_1= 0 
        i_2=0#erase!!!!!
        SFR_array = BC03dict['SFR']
        import random
        mock_gal_catalog=[]
        ##########################################erase!!!!!

				##AFter 97

                ##########################################erase!!!!!
                i_1=i_1+1  
                if i_1 in random_array:
                    i_2=i_2+1
                    zm=random.randint(1,20)*0.1
                    GAm = random.randint(2,4)
                    dm = z2Dlum(zm)
                    solarlum = const.L_sun.to(u.erg/u.second) #3.839e33
                    Nm = 10**GAm * 4* pi* dm**2 / (solarlum.value)/ (1+zm)
                    Nm =renorm_template('GA', Nm)
                    xm= gal_nu
                    ym= gal_Fnu_red * Nm #*solarlum.value *1e-23 *1e-3
                    SFRm = BC03dict['SFR'][:,agei,taui,:,:][0][0][0].value * Nm
                    #dlum= z2Dlum(1)
                    #lumfactor = (4. * pi * dlum**2.)

                    ###Choose data to simulate

                    ### look for this frequencies
                    #GALEX_1500|1516|15.296
                    #u_SDSS|3591|14.921
                    #g_SDSS|4723|14.802|
                    #r_SDSS|6213|14.683|
                    #i_SDSS|7523|14.600|
                    #z_SDSS|8855|14.529|
                    #F105W|10517|14.455
                    #F140W|13877|14.334
                    #F160W|15348|14.291
                    #H_2mass|16467|14.260
                    #Ks_2mass|21641|14.141
                    #IRAC1|35375|13.928|
                    #IRAC2|44751|13.826|

                    ### Calculate fluxes to extract from model. Redshift model and extract at telescope bands.
                    freq_obs_array=np.array([13.836, 13.928, 14.141, 14.260, 14.291, 14.331, 14.455, 14.529, 14.600, 14.683, 14.802, 14.921, 15.296])
                    freq_intrinsic_array= np.log10(10**freq_obs_array*(1+zm))
                    fidxs=[]
                    for fi in freq_intrinsic_array:        
                        fidx = (np.abs(np.log10(xm)-fi)).argmin()
                        fidxs.append(fidx)
                    fidxsa= np.array(fidxs)[0:2].tolist()
                    fidxsb= np.array(fidxs)[2:].tolist()
                    ymerr=np.hstack((ym[fidxsa]*0.5, ym[fidxsb]*0.1))

                    ### Construct output data, ID, SFR, z, nu, Fnu, Fnu_err
                    data =np.array([np.hstack((i_2,SFRm, zm, 10**freq_obs_array, ym[fidxs], ym[fidxs]*0.05))])

                    ### Do plots to compare
                    fig = plt.figure()
                    plt.errorbar(freq_obs_array, ym[fidxs]*(10**freq_obs_array), yerr= ymerr*(10**freq_obs_array), ecolor='k', linestyle='', ms=3)
                    plt.plot(np.log10(xm/(1+zm)),ym*(xm/(1+zm)),'r-',lw=0.5, label=' age='+str(np.log10(age_array.value[agei])) +\
                                                                                    ' tau='+str(tau_array.value[taui])+ \
                                                                                    ' EB-V='+str(ebvgal_array[ebvi]) +\
                                                                                    ' SFR='+str(SFRm)+', z='+str(zm))
                    plt.yscale('log')
                    plt.legend(prop={'size':7})
                    #plt.savefig('/Users/gcalistr/Documents/AGNfitter/data/mockGalSED2_'+str(i_2)+'.png')
                    plt.close

                    mock_gal_catalog.append(data)
                
        mock_gal_catalog= np.reshape(np.array(mock_gal_catalog),[42,42])

        ### Write mock catalog
        # ascii.write(mock_gal_catalog, '/Users/gcalistr/Documents/AGNfitter/data/mockGalSED2.dat', names=['N','SFR', 'z', \
        #                                                                                                 'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12','x13',\
        #                                                                                                 'y1','y2','y3','y4','y5','y6','y7','y8','y9','y10','y11','y12','y13',\
        #                                                                                                 'ye1','ye2','ye3','ye4','ye5','ye6','ye7','ye8','ye9','ye10','ye11','ye12','ye13'],overwrite=True)
        #########################################

                    