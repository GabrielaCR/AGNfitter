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

    cat['path'] = '/home/user/AGNfitter/' #path to the AGNfitter code
    cat['filename'] = cat['path']+ 'data/catalog_example.txt'  
    cat['filetype'] = 'ASCII'		## catalog file type: 'ASCII' or 'FITS'. 
    cat['name'] = 0			## If ASCII: Column index (int) of source IDs
                                   	## If FITS : Column name (str). E.g. 'ID'
    cat['redshift'] = 1       		## If ASCII:  Column index(int) of redshift 
                                    	## If FITS : Column name (str). E.g. z'

   ##FREQUENCIES/WAVELENGTHS 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['freq/wl_list'] = np.arange(2,147,3).tolist()                         

                                        ## If ASCII: List of column indexes (int), 

                                        ##           corresponding to freq/wl.                                  
    #cat['freq/wl_suffix'] = '_wl'      ## If FITS: common ending to wavelength column names

    cat['use_central_wavelength'] = True # Option to use central wavelength if no wavelengths in table

    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*
                                         ## 'frequency' or 'wavelength'?
    cat['freq/wl_unit'] = u.micron	 ## Astropy unit of freq or wavelength


    ##FLUXES 
    ## if ASCII specify 'freq/wl_list', if FITS specify 'freq/wl_suffix'
    cat['flux_in_magAB'] = False 	# Option to calculate flux and flux_error from magnitude AB.
    cat['flux_unit'] = u.Jy * 1e-3	## Astropy unit of *flux* (astropy-units)
    cat['flux_list'] = np.arange(3,148,3).tolist()     
                                        ## If ASCII: List of column indexes (int)
    #cat['flux_suffix'] = '_f'          ## If FITS: Common ending of all flux column names (str)    
    cat['fluxerr_list'] = np.arange(4,149,3).tolist()
                                        ## If ASCII: List of column indexes (int)
    #cat['fluxerr_suffix'] = '_e'       ## If FITS: common ending to fluxerr column names (str)
    cat['err+10%flux_moreflex'] = False
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

    filters['dict_zarray'] = np.array([]) # Deprecated. The grid of redshifts needed to fit your catalog
    filters['path'] = 'models/FILTERS/' 
    filters['filterset'] = 'example_46datapointa' ## 'filterset_default' (for the test case),
                                               ## for the user's case: customize, eg. filtersv1

    filters['4XMM_EP5'] = [True, 48]  #4XMM_EP5
    filters['4XMM_EP4'] = [True, 47]  #4XMM_EP4
    filters['4XMM_EP3'] = [True, 46]  #4XMM_EP3
    filters['4XMM_EP2'] = [True, 45]  #4XMM_EP2
    filters['RAD_1_4GHz'] = [True, 44]   #1_4GHz
    filters['RAD_5GHz'] = [True, 43]    #5GHz
    filters['RAD_15GHz'] = [True, 42]  #RAD_15GHz 
    filters['SPIRE500'] = [True, 41]   #S500
    filters['SPIRE350'] = [True, 40]   #S350
    filters['SPIRE250'] = [True, 39]   #S250
    filters['PACS160'] = [True, 38]    #P160
    filters['PACS100'] = [True, 37]    #P100
    filters['PACS70'] = [True, 36]     #P70
    filters['M1'] = [True, 35]         #M1
    filters['WISE4'] = [True, 34]      #W4
    filters['WISE3'] = [True, 33]      #W3
    filters['IRAC4'] = [True, 32]      #I4
    filters['IRAC3'] = [True, 31]      #I3
    filters['WISE2'] = [True, 30]      #W2
    filters['IRAC2'] = [True, 29]      #I2
    filters['IRAC1'] = [True, 28]      #I1
    filters['WISE1'] = [True, 27]      #W1
    filters['Ks_2mass'] = [True, 26]   #Ks
    filters['H_2mass'] = [True, 25]    #H
    filters['J_2mass'] = [True, 24]    #J
    filters['y_PS1_good'] = [True, 23] #py  y_PS1_good
    filters['SKY_z'] = [True, 22]      #Skz
    filters['z'] = [True, 21]          #z
    filters['z_PS1'] = [True, 20]      #pz
    filters['SKY_i'] = [True, 19]      #Ski
    filters['i_PS1'] = [True, 18]      #pi
    filters['i'] = [True, 17]          #i
    filters['SKY_r'] = [True, 16]      #Skr
    filters['r_PS1'] = [True, 15]      #pr
    filters['r'] = [True, 14]          #r
    filters['SWIFT_V'] = [True, 13]    #SV
    filters['SKY_g'] = [True, 12]      #Skg
    filters['g_PS1'] = [True, 11]      #pg
    filters['g'] = [True, 10]          #g
    filters['SWIFT_B'] = [True, 9]     #SB
    filters['SKY_v'] = [True, 8]       #Skv
    filters['u'] = [True, 7]           #u
    filters['SKY_u'] = [True, 6]       #Sku
    filters['SWIFT_U'] = [True, 5]     #SU
    filters['SWIFT_W1'] = [True, 4]    #SW1
    filters['GALEX_2500'] = [True, 3]  #NUV
    filters['SWIFT_M2'] = [True, 2]    #SM2
    filters['SWIFT_W2'] = [True, 1]    #SW2
    filters['GALEX_1500'] = [True, 0]  #FUV


    filters['add_filters']= False # If 'True' please add them below in ADD FILTERS

    """==================================
    ADD FILTERS (optional)
    =================================="""

    ADDfilters=dict()
    ADDfilters['names'] = []		## (string/list of strings)User especified filter names. 
                                	## If name has been previously used, an error message will appear. 
    ADDfilters['filenames'] = []	## (string/list of strings) File names of the new filters. 
                                	## File format: 2 columns of 1) freq/wavelength 2) Throughput. 
                                	## Path assumed is the cat['path'] especified above. 
                                	## Example: 'models/FILTERS/my_new_filter.txt'

    ADDfilters['freq/wl_format'] = ['wavelength'] * len(ADDfilters['names']) 	## Info about the column 1 of your filter file.
                                                                             	## Options: 'wavelength' or 'frequency'.    
    ADDfilters['freq/wl_unit'] =  [u.Angstrom]* len(ADDfilters['names']) 	## (Astropy Unit) Info about the column 1 
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


    models['GALAXY'] = 'BC03_metal'   	### Current options:
                                	### 'BC03' (Bruzual & Charlot 2003)
                                	### 'BC03_metal' (Bruzual & Charlot 2003), with metallicities
    models['STARBURST'] = 'S17' 	### Current options:
                                	### 'DH02_CE01' (Dale & Helou 2002 + Chary & Elbaz 2001)
                                	### 'S17' (Schreiber et al. 2017 (submitted)) with very small grains (VSG) correction 
                                	### 'S17_radio' Previous model + the IR-radio correlation (Bell,2003) for the host galaxy radio emission

    models['BBB'] = 'THB21' 	### Current options:
                         	### 'R06' (Richards et al. 2006) 
                         	### 'SN12' (Slone&Netzer 2012)
                         	### 'KD18' (Kubota and Done 2018)
                         	### 'THB21' (Temple et al. 2021)

    models['TORUS'] ='SKIRTOR_1P' ### Current options:
                           	### 'S04' (Silva et al. 2004)
                           	### 'NK0' (Nenkova et al. 2008)
                           	### 'NK0_2P' (Nenkova et al. 2008) with averaged SEDs for each inclination and openning angle
                           	### 'NK0_3P' (Nenkova et al. 2008) with averaged SEDs for each inclination, openning angle and optical depth
                           	### 'SKIRTOR_1P' (Stalevski et al. 2016) model with averaged SEDs for each inclination
                           	### 'SKIRTOR_2P' (Stalevski et al. 2016) model with averaged SEDs for each inclination and openning angle
                           	### 'SKIRTOR_3P' (Stalevski et al. 2016) model with averaged SEDs for each inclination, openning angle and optical depth
                           	### 'CAT3D_3P' CAT3D-Wind model with averaged SEDs for each inclination angle, a and fwd parameters.

    models['XRAYS'] = False 	### If X-ray data is available and informative for the fit

    models['RADIO'] = False 	### If radio data is available and informative for the fit. The code choose the radio model according 
				### to the number of available bands.
                                ### Simple power law with slope= -0.75 for 1 band or SPL with slope as a fitting parameter for 2 bands 
                                ### Double power law for 3 bands or more. The fitting parameters are alpha1-alpha2, nu_t for 3 bands or 
				### alpha1, alpha2 and nu_t for 4 bands.

    models['Blazar'] = 'None'    ### Same previous models but now cutoff energy is included as an extra fitting parameters.
				### Options: 'SPL', 'DPL-3', 'DPL-4' for a simple power law and double power law with 3 and 4 fitting parameters
                                ### 'None' for non-blazar fit

    models['PRIOR_energy_balance'] = 'Flexible' ### Default:'Flexible'
                                          	### 'Flexible': Sets a lower limit to the dust emission luminosity ('starburst' model)
                                          	### as given by the observed attenuation in the stellar component SED.
                                          	### 'Restrictive': Favours a combination of parameters such that the luminosity of the cold dust 
                                          	### and that attenuated in the stellar component are equal.

    models['PRIOR_AGNfraction'] = True  ### Default: True
                                        ### True: - *IF* blue/UV bands (around 1500 Angstrom) are 10 times higher than expected by the galaxy
                                        ### luminosity function by Parsa, Dunlop et al. 2014, AGN-to-GAL ratios lower than 1 (log =0) are rejected. 
					### It then applies a Gaussian prior probability with log ratio=2, with a sigma of 2.
                                        ### - In this cases it also applies a Gaussian prior on the galaxy normalization, i.e. stellar mass (usually 
					###   unconstrained in these cases) to populate physically expected ranges for QSO hosts -> 10^9 - 10^11. 
                                        ### - *ELSE IF* blue/UV bands (around 1500 Angstrom) are below 10 times the expected value by Parsa, Dunlop et 
					### al. 2014. this option gives preference to galaxy contribution in the optical UV, with Gaussian prior 
					###  probability centered on AGN to GALAXY log ratios of -1. and sigma 1, i.e. accretion disk is disfavoured at least the data strongly prefers it.
                                        ### False:- Non-informative prior

    models['PRIOR_midIR_UV'] = False
    models['PRIOR_galaxy_only'] = False ### Default:False 
                                        ### True: sets all AGN contribution to 0.ß

    models['turn_on_AGN'] = True	### Includes accretion disk, torus, X-rays and radio components is AGN is turn on. 
					### Otherwise just fit stellar and cold dust emission.

    return models

def MCMC_settings():

    """==================================
    Set your preferences for the MCMC sampling.
    =================================="""

    mc = dict()
    mc['sampling_algorithm'] = 'ultranest' #ultranest or emcee

    # If mcmc algorithm is emcee, please define the following values, otherwise just keep the current values
    mc['Nwalkers'] = 100	## number of walkers 
    mc['Nburnsets']= 2		## number of burn-in sets
    mc['Nburn'] = 12000		## length of each burn-in sets 
    mc['Nmcmc'] = 25000		## length of each burn-in sets 
    mc['iprint'] = 1000		## show progress in terminal in steps of this many samples

   # If mcmc algorithm is ultranest, please define the following values, otherwise just keep the current values
    mc['direction_generation'] = 'de-mix' # Options: 1) 'de-mix' for mixture random direction
                                          #          2) 'region-slice' for region oriented vector
                                          #          3) 'cube-ortho-harm' for random orthogonal direction vector
    mc['live_points'] = 400	# Number of initial live points to explore the space parameter
    mc['min_ess'] = 400		# Minimum number of effective samples
    mc['num_loops'] = 0		# How many times to go back and improve

    return mc

def OUTPUT_settings():

    """==================================
    Set your preferences for the production of OUTPUT files. 
    =================================="""

    out = dict()

    out['plot_format'] = 'pdf'            #pdf
    out['band_indicators'] = True         #Show tags of the wavelength (e.g.'Radio', 'FIR')
    out['realizations_format'] = 'shaded' #'shaded' for shaded area given by maximum and minimum fluxes of 100 random realizations 
                                          # or 'curves' to plot 10 random model combinations (like Calistro Rivera 2016)

    #CHAIN TRACES
    out['plot_tracesburn-in'] = False    
    out['plot_tracesmcmc'] = True

    #BASIC OUTPUT
    out['Nsample'] = 1000 ## out['Nsample'] * out['Nthinning'] <= out['Nmcmc']
    out['Nthinning'] = 10 ## This describes thinning of the chain to sample
    out['writepar_meanwitherrors'] = True ##Write output values for all parameters in a file.
    out['plot_posteriortriangle'] = True ##Plot triangle with all parameters' PDFs?

    #INTEGRATED LUMINOSITIES
    out['calc_intlum'] = True  
    out['save_posterior_luminosities']= False
    out['save_posteriors']= True
    out['realizations2int'] = 100 #This process is very time consuming.
                                #Around 100-1000 is recomendend for computational reasons.
                                #If you want to plot posterior triangles of 
                                #the integrated luminosities, should be > 1000.
    out['plot_posteriortrianglewithluminosities'] = False  # requires out['calc_intlum']=True 

    #INTEGRATION RANGES

    out['intlum_models'] = ['sb','bbb', 'bbbdered', 'gal', 'AGNfrac_opt', 'tor','tor+bbb','AGNfrac_IR', 'gal', 'bbb', 'AGNfrac_opt','tor','sb']  
	#You can add for radio 'AGNfrac_rad', 'Lx0.5-2keV','Lx2-8keV' with wavelenghts in microns 

    out['intlum_freqranges_unit'] = u.micron   #Astropy unit 

    out['intlum_freqranges'] = np.array([[8.,1000.],[0.1,1.],[0.1,1.],[0.1,1.], [0.1,1.], [1.,30.], [0.1, 30], [8.,1000.],[0.4, 0.5],[0.4, 0.5], [0.4, 0.5], [6,6], [1.,30.]])

    out['intlum_names'] = ['LIR(8-1000)','Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'AGNfrac(0.1-1)', 'Ltor(1-30)', 'LAGN(0.1-30)', 'AGNfrac(8-1000)', 'Lgal(0.4-0.5)', 'Lbbb(0.4-0.5)', 'AGNfrac(0.4-0.5)','Ltor(6)', 'Lsb(1-30)']


    #SED PLOTTING
    out['realizations2plot'] = 10 #For plotting curves we suggest 10 and for shaded areas 100
    out['plot_best_fit'] = True
    out['plot_median_fit'] = False 
    out['plot_residuals']= True

    out['saveSEDresiduals'] = True
    out['plotSEDrealizations'] = True
    out['saveSEDrealizations'] = True

    return out
