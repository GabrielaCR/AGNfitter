
"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DICTIONARIES_AGNFitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains all functions which are needed to construct the total model of AGN. 

##For constructing a new dictionary, 
(in cases: 1)add a filter which is not included, 
2) need finer grid for better S/N data)
 see DICTIONARIES_AGNfitter.py  
    

"""
import numpy as np
import sys
from collections import defaultdict

import MODEL_AGNfitter as model
import time
import cPickle
import shelve
from astropy import units as u 


class MODELSDICT:


	"""
	Class MODELSDICT

	Builds a dictionary of model templates. 

	##input: 
	- filename of the dictionary you want to create
	- the path whre it will be located

	- Also variables self.ebvgal_array,self.ebvbbb_array, self.z_array
	  can be change by the user, for a finer grid in this parameters.
	##bugs: 

	"""     

	def __init__(self, filename, path):
		self.filename = filename
		self.path=path
		self.ebvgal_array = np.array(np.arange(0.,100.,5.)/100)
		self.ebvbbb_array = np.array(np.arange(0.,100.,5.)/100)
		self.z_array =np.array([0.14,1.58])

	def build(self):

		f = open(self.filename, 'wb')

		COSMOS_modelsdict = dict()

		print 'MODELSDICT.build'
		print 'Constructing Dictionary of models.' 
		print '--------------------------------------'
		print 'Make sure the filterset contains all the photometric bands'
		print 'needed by your catalog.'
		print 'This process might take a while, but you have to do it only once.'
		print 'If you interrupt it, please trash the empty file created.'
		print ''

		for z in self.z_array:#use this range to make the 'z_array_in_dict' in RUN_AGNfitter.   				
			filterset='BANDSET_default'	# See function filter_dictionaries()	
			filterdict = filter_dictionaries(filterset, self.path)
			dict_modelsfiltered = self.construct_dictionaryarray_filtered(z, filterdict, self.path)
			COSMOS_modelsdict[str(z)] = dict_modelsfiltered
   	
		print 'Dictionary has been created in :', self.filename

		cPickle.dump(COSMOS_modelsdict, f, protocol=2)
		f.close()



	def construct_dictionaryarray_filtered(self, z, filterdict,path):

		"""
		Construct the dictionaries of fluxes at bands (to compare to data), 
		and dictionaries of fluxes over the whole spectrum, for plotting.
		"""

		GALAXYFdict_filtered = dict()
		GALAXY_SFRdict = dict()
		STARBURSTFdict_filtered = dict()		
		BBBFdict_filtered = dict()
		TORUSFdict_filtered = dict()

		GALAXYFdict_4plot = dict()
		STARBURSTFdict_4plot = dict()		
		BBBFdict_4plot = dict()
		TORUSFdict_4plot = dict()



		#OPENING TEMPLATES AND BUILDING DICTIONARIES

		#Call object containing all galaxy models 	
		galaxy_object = cPickle.load(file(path + 'models/GALAXY/bc03_275templates.pickle', 'rb')) 
		_, ageidx, tauidx, _, _,_ =  np.shape(galaxy_object.SED)
		#Construct dictionaries 
		for taui in range(tauidx):
			for agei in range(ageidx):

				gal_wl, gal_Fwl =  galaxy_object.wave, galaxy_object.SED[:,agei,taui,:,:,:].squeeze()
				gal_nus= gal_wl.to(u.Hz, equivalencies=u.spectral())[::-1]#invert
				gal_Fnu= (gal_Fwl * 3.34e-19 * gal_wl**2.)[::-1]  
				gal_SFR= galaxy_object.SFR[:,agei,taui,:,:].squeeze()
				GALAXY_SFRdict[str(galaxy_object.tau.value[taui]),str(galaxy_object.tg.value[agei])] = gal_SFR

				for EBV_gal in self.ebvgal_array:
					#Apply reddening				
					gal_nu, gal_Fnu_red = model.GALAXY_nf2(gal_nus.value[0:len(gal_nus):3], gal_Fnu.value[0:len(gal_nus):3], EBV_gal)	
					GALAXYFdict_4plot[str(galaxy_object.tau.value[taui]),str(galaxy_object.tg.value[agei]), str(EBV_gal)] = \
																								np.log10(gal_nu), gal_Fnu_red
					#Projection of filter curves on models
					bands,  gal_Fnu_filtered =  model.filters1(np.log10(gal_nu), gal_Fnu_red, filterdict, z)			
					GALAXYFdict_filtered[str(galaxy_object.tau.value[taui]),str(galaxy_object.tg.value[agei]), str(EBV_gal)] = \
																										bands, gal_Fnu_filtered



		#Call object containing all starburst models 	
		starburst_object = cPickle.load(file(path + 'models/STARBURST/dalehelou_charyelbaz_v1.pickle', 'rb')) 
		irlumidx = len(starburst_object.SED)
		#Construct dictionaries 
		for irlumi in range(irlumidx):
			sb_nu0, sb_Fnu0 = starburst_object.wave[irlumi], starburst_object.SED[irlumi].squeeze()
			STARBURSTFdict_4plot[str(starburst_object.irlum[irlumi])] = sb_nu0, sb_Fnu0
			bands, sb_Fnu_filtered = model.filters1(sb_nu0, sb_Fnu0, filterdict, z)
			STARBURSTFdict_filtered[str(starburst_object.irlum[irlumi])] = bands, sb_Fnu_filtered
			if np.amax(sb_Fnu_filtered) == 0:
				print 'Error: something is wrong in the calculation of STARBURST flux'



		#No object to call since bbb is only one model 	
		bbb_object = cPickle.load(file(path + 'models/BBB/richards.pickle', 'rb')) 

		bbb_nu, bbb_Fnu = bbb_object.wave, bbb_object.SED.squeeze()
		#Construct dictionaries
		for EBV_bbb in self.ebvbbb_array:
			bbb_nu0, bbb_Fnu_red = model.BBB_nf2(bbb_nu, bbb_Fnu, EBV_bbb, z )
			BBBFdict_4plot[str(EBV_bbb)] =bbb_nu0, bbb_Fnu_red
			bands, bbb_Fnu_filtered = model.filters1(bbb_nu0, bbb_Fnu_red, filterdict,z)
			BBBFdict_filtered[str(EBV_bbb)] = bands, bbb_Fnu_filtered
			if np.amax(bbb_Fnu_filtered) == 0:
				print 'Error: something is wrong in the calculation of BBB flux'			



		#Call object containing all torus models 	
		torus_object = cPickle.load(file(path + 'models/TORUS/silva_v1.pickle', 'rb')) 
		nhidx=len(torus_object.SED)
		#Construct dictionaries 
		for nhi in range(nhidx):

			tor_nu0, tor_Fnu0 = torus_object.wave[nhi], torus_object.SED[nhi].squeeze()
			TORUSFdict_4plot[str(torus_object.nh[nhi])] = tor_nu0, tor_Fnu0

			bands, tor_Fnu_filtered = model.filters1(tor_nu0, tor_Fnu0, filterdict, z)
			TORUSFdict_filtered[str(torus_object.nh[nhi])] = bands, tor_Fnu_filtered
			if np.amax(tor_Fnu_filtered) == 0:
				print 'Error: something is wrong in the calculation of TORUS flux'




		return STARBURSTFdict_filtered , BBBFdict_filtered, GALAXYFdict_filtered, TORUSFdict_filtered, \
			   STARBURSTFdict_4plot , BBBFdict_4plot, GALAXYFdict_4plot, TORUSFdict_4plot,GALAXY_SFRdict
			   




def dictkey_arrays(MODELSdict):

	"""
	Construct the dictionaries of fluxes at bands (to campare to data), 
	and dictionaries of fluxes over the whole spectrum, for plotting.

	##input:

	##output:
	"""

	STARBURSTFdict , BBBFdict, GALAXYFdict, TORUSFdict, _,_,_,_,GALAXY_SFRdict= MODELSdict
	tau_dict= np.array(list(GALAXYFdict.keys()))[:,0]
	age_dict= np.array(list(GALAXYFdict.keys()))[:,1]
	ebvg_dict = np.array(list(GALAXYFdict.keys()))[:,2]

	irlum_dict = np.array(list(STARBURSTFdict.keys()))
	nh_dict = np.array(list(TORUSFdict.keys()))
	ebvb_dict = np.array(list(BBBFdict.keys()))


	#For computational reasons (to be used in PARAMETERspace_AGNfitter.py)
	class gal_class:
		def __init__(self, tau_dict, age_dict, ebvg_dict):
			self.tau_dict =tau_dict
			self.age_dict= age_dict
			self.ebvg_dict = ebvg_dict
			self.tau_dict_float =tau_dict.astype(float)
			self.age_dict_float= age_dict.astype(float)
			self.ebvg_dict_float = ebvg_dict.astype(float)

		def nearest_par2dict(self, tau, age, ebvg):	
			taui =np.abs(self.tau_dict_float-tau).argmin()
			agei= np.abs(self.age_dict_float-age).argmin()
			ebvgi = np.abs(self.ebvg_dict_float-ebvg).argmin()
			self.t = tau_dict[taui]
			self.a= age_dict[agei]
			self.e= ebvg_dict[ebvgi]



	gal_obj = gal_class(tau_dict, age_dict, ebvg_dict)

	return gal_obj, irlum_dict, nh_dict, ebvb_dict, GALAXY_SFRdict





def filter_dictionaries(filterset, path):

	"""
	Constructs the dictionaries of fluxes 
	1) specifically for your photometric bands (to campare to data), and
	2) dictionaries of fluxes for the whole spectrum, for plotting.

	input
	-------
	- filterset: Here we have two types of filterset: 
	'COSMOS1' or 'BANDSET_default',
	which comprise together most of the existent bands.

	If the bands you use are inside this sets, 
	but you don't have all, you can still use them, no need to change.

	If some band is missing, you can construct an own filterdictionary
	with the CLASS MODELSDICT. Create your own filterset,
	combining the filters we have, and adding your own ones.
	**The bands must have the  order correspondingto the elements in the lists**
	Change the name of the filterset in 
	RUN_AGNfitter_multi.py, 
	function catalog_settings, cat['dict_path'] = 'myfilterset'.

	dependency
	----------
	This function is called in the CLASS MODELSDICT

	"""

	if filterset == 'COSMOS1':

		bands = [ 12.27250142,  12.47650142,  12.6315,  13.09650142,  13.59050142, 13.72250142,  13.82750142,  13.92950142,  14.14450142, 14.26450142, 14.38150142 , 14.52150142 , 14.59450142  ,14.68250142,  14.74050142, 14.80250142  ,14.82950142 , 14.88450142  ,15.11350142, 15.28650142]

		#INFRAROJO

		H160band_file = path + 'models/FILTERS/HERSCHEL/PACS_160mu.txt'
		H160_lambda, H160_factor =  np.loadtxt(H160band_file, usecols=(0,1),unpack= True)

		H100band_file =path + 'models/FILTERS/HERSCHEL/PACS_100mu.txt'
		H100_lambda, H100_factor =  np.loadtxt(H100band_file, usecols=(0,1),unpack= True)

		M70band_file = path + 'models/FILTERS/SPITZER/mips70.res'
		M70_lambda, M70_factor =  np.loadtxt(M70band_file, usecols=(0,1),unpack= True)

		M24band_file =  path + 'models/FILTERS/SPITZER/mips24.res'
		M24_lambda, M24_factor =  np.loadtxt(M24band_file, usecols=(0,1),unpack= True)

		I4band_file = path + 'models/FILTERS/SPITZER/irac_ch4.res'
		I4_lambda, I4_factor =  np.loadtxt(I4band_file, usecols=(0,1),unpack= True)

		I3band_file =  path + 'models/FILTERS/SPITZER/irac_ch3.res'
		I3_lambda, I3_factor =  np.loadtxt(I3band_file, usecols=(0,1),unpack= True)

		I2band_file = path + 'models/FILTERS/SPITZER/irac_ch2.res'
		I2_lambda, I2_factor =  np.loadtxt(I2band_file, usecols=(0,1),unpack= True)

		I1band_file = path + 'models/FILTERS/SPITZER/irac_ch1.res'
		I1_lambda, I1_factor =  np.loadtxt(I1band_file, usecols=(0,1),unpack= True)

		Kband_file = path + 'models/FILTERS/2MASS/Ks_2mass.res'
		K_lambda, K_factor =  np.loadtxt(Kband_file, usecols=(0,1),unpack= True)

		Hband_file = path + 'models/FILTERS/2MASS/H_2mass.res'
		H_lambda, H_factor =  np.loadtxt(Hband_file, usecols=(0,1),unpack= True)

		Jband_file = path + 'models/FILTERS/2MASS/J_2mass.res'
		J_lambda, J_factor =  np.loadtxt(Jband_file, usecols=(0,1),unpack= True)

		zband_file =path + 'models/FILTERS/SUBARU/z_subaru.res'
		z_lambda, z_factor =  np.loadtxt(zband_file, usecols=(0,1),unpack= True)

		iband_file = path + 'models/FILTERS/CHFT/i_megaprime_sagem.res'
		i_lambda, i_factor =  np.loadtxt(iband_file, usecols=(0,1),unpack= True)

		rband_file = path + 'models/FILTERS/SUBARU/r_subaru.res'
		r_lambda,r_factor =  np.loadtxt(rband_file, usecols=(0,1),unpack= True)

		Vband_file = path + 'models/FILTERS/SUBARU/V_subaru.res'
		V_lambda, V_factor =  np.loadtxt(Vband_file, usecols=(0,1),unpack= True)

		gband_file =path + 'models/FILTERS/SUBARU/g_subaru.res'
		g_lambda,g_factor =  np.loadtxt(gband_file, usecols=(0,1),unpack= True)

		Bband_file = path + 'models/FILTERS/SUBARU/B_subaru.res'
		B_lambda, B_factor =  np.loadtxt(Bband_file, usecols=(0,1),unpack= True)

		uband_file = path + 'models/FILTERS/CHFT/u_megaprime_sagem.res'
		u_lambda, u_factor =  np.loadtxt(uband_file, usecols=(0,1),unpack= True)

		NUVband_file = path + 'models/FILTERS/GALEX/galex2500.res'
		NUV_lambda, NUV_factor =  np.loadtxt(NUVband_file, usecols=(0,1),unpack= True)

		FUVband_file = path + 'models/FILTERS/GALEX/galex1500.res'
		FUV_lambda, FUV_factor =  np.loadtxt(FUVband_file, usecols=(0,1),unpack= True)

		
		files = [H160band_file, H100band_file, M70band_file, M24band_file,I4band_file , \
				I3band_file, I2band_file, I1band_file, Kband_file, Hband_file, Jband_file,\
				zband_file  , iband_file, rband_file, Vband_file, gband_file , Bband_file,\
				uband_file, NUVband_file, FUVband_file]

		lambdas = [H160_lambda, H100_lambda, M70_lambda, M24_lambda, I4_lambda, I3_lambda, \
				I2_lambda, I1_lambda,K_lambda,H_lambda , J_lambda,z_lambda, i_lambda,r_lambda,\
				V_lambda ,g_lambda , B_lambda, u_lambda, NUV_lambda,FUV_lambda]

		factors = [H160_factor, H100_factor, M70_factor, M24_factor, I4_factor, I3_factor,\
				 I2_factor, I1_factor, K_factor, H_factor , J_factor, z_factor, i_factor,\
				  r_factor ,V_factor ,g_factor , B_factor , u_factor , NUV_factor ,FUV_factor ]

        #dictionaries lambdas_dict, factors_dict

		files_dict = defaultdict(list)
		lambdas_dict = defaultdict(list)
		factors_dict = defaultdict(list)

		for i in range(len(files)):

			files_dict[bands[i]].append(files[i])
			lambdas_dict[bands[i]].append(lambdas[i])
			factors_dict[bands[i]].append(factors[i])

				
	if filterset == 'BANDSET_default':

		#central frequencies in log10(nu)
		bands = [ 11.77815 , 11.933053, 12.0791812, 13.09650142,  13.59050142, 13.69897, 13.80919, 13.90609, 14.12996,  14.2499 , 14.3662, 14.4491, 14.5169, 14.5748, 14.6686,  14.8239,  14.91335, 15.1135]


		#All filter files
		H500band_file = path + 'models/FILTERS/HERSCHEL/SPIRE_500mu.txt'
		H500_lambda, H500_factor =  np.loadtxt(H500band_file, usecols=(0,1),unpack= True)

		H350band_file = path + 'models/FILTERS/HERSCHEL/SPIRE_350mu.txt'
		H350_lambda, H350_factor =  np.loadtxt(H350band_file, usecols=(0,1),unpack= True)

		H250band_file = path + 'models/FILTERS/HERSCHEL/SPIRE_250mu.txt'
		H250_lambda, H250_factor =  np.loadtxt(H250band_file, usecols=(0,1),unpack= True)

		M24band_file =  path + 'models/FILTERS/SPITZER/mips24.res'
		M24_lambda, M24_factor =  np.loadtxt(M24band_file, usecols=(0,1),unpack= True)


		I4band_file = path + 'models/FILTERS/SPITZER/irac_ch4.res'
		I4_lambda, I4_factor =  np.loadtxt(I4band_file, usecols=(0,1),unpack= True)

		I3band_file =  path + 'models/FILTERS/SPITZER/irac_ch3.res'
		I3_lambda, I3_factor =  np.loadtxt(I3band_file, usecols=(0,1),unpack= True)

		I2band_file = path + 'models/FILTERS/SPITZER/irac_ch2.res'
		I2_lambda, I2_factor =  np.loadtxt(I2band_file, usecols=(0,1),unpack= True)

		I1band_file = path + 'models/FILTERS/SPITZER/irac_ch1.res'
		I1_lambda, I1_factor =  np.loadtxt(I1band_file, usecols=(0,1),unpack= True)

		
		Kband_file = path + 'models/FILTERS/2MASS/Ks_2mass.res'
		K_lambda, K_factor =  np.loadtxt(Kband_file, usecols=(0,1),unpack= True)

		Hband_file = path + 'models/FILTERS/2MASS/H_2mass.res'
		H_lambda, H_factor =  np.loadtxt(Hband_file, usecols=(0,1),unpack= True)

		Jband_file = path + 'models/FILTERS/2MASS/J_2mass.res'
		J_lambda, J_factor =  np.loadtxt(Jband_file, usecols=(0,1),unpack= True)

		Yband_file = path + 'models/FILTERS/VISTA/Y_uv.res'
		Y_lambda, Y_factor =  np.loadtxt(Yband_file, usecols=(0,1),unpack= True)

		zband_file =path + 'models/FILTERS/SUBARU/z_subaru.res'
		z_lambda, z_factor =  np.loadtxt(zband_file, usecols=(0,1),unpack= True)

		iband_file = path + 'models/FILTERS/CHFT/i_megaprime_sagem.res'
		i_lambda, i_factor =  np.loadtxt(iband_file, usecols=(0,1),unpack= True)

		rband_file = path + 'models/FILTERS/SUBARU/r_subaru.res'
		r_lambda,r_factor =  np.loadtxt(rband_file, usecols=(0,1),unpack= True)
		
		Bband_file = path + 'models/FILTERS/SUBARU/B_subaru.res'
		B_lambda, B_factor =  np.loadtxt(Bband_file, usecols=(0,1),unpack= True)

		uband_file = path + 'models/FILTERS/CHFT/u_megaprime_sagem.res'
		u_lambda, u_factor =  np.loadtxt(uband_file, usecols=(0,1),unpack= True)

		NUVband_file = path + 'models/FILTERS/GALEX/galex2500.res'
		NUV_lambda, NUV_factor =  np.loadtxt(NUVband_file, usecols=(0,1),unpack= True)

		#List of file names
		files = [ H500band_file, H350band_file, H250band_file, M24band_file, I4band_file ,\
		 I3band_file, I2band_file, I1band_file, Kband_file, Hband_file, Jband_file,\
		  Yband_file, zband_file , iband_file, rband_file,  Bband_file,  uband_file, \
		  NUVband_file]

		#List of all lambdas
		lambdas = [H500_lambda, H350_lambda, H250_lambda, M24_lambda, I4_lambda , I3_lambda, \
		I2_lambda, I1_lambda,  K_lambda, H_lambda, J_lambda, Y_lambda,  z_lambda, i_lambda, \
		r_lambda, B_lambda,  u_lambda, NUV_lambda]

		#list of all factors corresponding to the lambdas
		factors = [ H500_factor, H350_factor, H250_factor, M24_factor, I4_factor , I3_factor, \
		I2_factor, I1_factor, K_factor, H_factor, J_factor, Y_factor, z_factor, i_factor, \
		r_factor,  B_factor,  u_factor, NUV_factor]


        #dictionaries lambdas_dict, factors_dict
		files_dict = defaultdict(list)
		lambdas_dict = defaultdict(list)
		factors_dict = defaultdict(list)

		for i in range(len(files)):

			files_dict[bands[i]].append(files[i])
			lambdas_dict[bands[i]].append(lambdas[i])
			factors_dict[bands[i]].append(factors[i])

		
	return bands, files_dict, lambdas_dict, factors_dict




