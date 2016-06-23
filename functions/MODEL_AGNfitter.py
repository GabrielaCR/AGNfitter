

"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             MODEL_AGNfitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains all functions which are needed to construct the total model of AGN. 
The functions here translate the parameter space points into total fluxes dependin on the models chosen.

Functions contained here are the following:

pick_STARBURST_template
pick_GALAXY_template
pick_TORUS_template
pick_EBV_grid


STARBURST_nf
BBB_nf
GALAXY_nf
TORUS_nf

"""

import numpy as np
from math import exp,pi, sqrt
import matplotlib.pyplot as plt
import time
from scipy.interpolate import interp1d
from scipy.integrate  import quad, trapz
import astropy.constants as const
import astropy.units as u




"""==============================
PICKING TEMPLATES
==============================

Functions which compensate for the discreteness of pur models. 
It infers the existent model dictionary key 'par_key',
from the continous valus par_mcmc, through NearestNeighbour interpolation.

"""


def pick_STARBURST_template(ir_lum, irlum_dict):

	idx = (np.abs(irlum_dict.astype(float)-ir_lum)).argmin()
	return irlum_dict[idx]

def pick_BBB_template(ebvb,ebvb_dict):
	
	ebvb_idx = (np.abs(ebvb_dict.astype(float)-ebvb)).argmin()
	return ebvb_dict[ebvb_idx]

def pick_GALAXY_template( tau, age, ebvg, tau_dict, age_dict, ebvg_dict):
	tauidx = (np.abs(tau_dict.astype(float)-tau)).argmin()	
	ageidx = (np.abs(age_dict.astype(float)-age)).argmin()
	ebvidx = (np.abs(ebvg_dict.astype(float)-ebvg)).argmin()
	return tau_dict[tauidx], age_dict[ageidx], ebvg_dict[ebvidx]

def pick_TORUS_template(nh, nh_dict):

	idx = (np.abs(nh_dict.astype(float)-nh)).argmin()
	return nh_dict[idx]

def pick_EBV_grid (EBV_array, EBV):

	idx = (np.abs(EBV_array-EBV)).argmin()
	EBV_fromgrid  = EBV_array[idx]

	return EBV_fromgrid

#==============================
# MAXIMAL POSSIBLE AGE FOR GALAXY MODEL
#==============================


def maximal_age(z):

	z = np.double(z)
	#Cosmological Constants	
	O_m = 0.266
	O_r =  0.
	O_k= 0.
	O_L = 1. - O_m
	H_0 = 74.3 #km/s/Mpc
	H_sec = H_0 / 3.0857e19 
	secondsinyear = 31556926
	ageoftheuniverse = 13.798e9

	# Equation for the time elapsed since z and now

	a = 1/(1+z)
	E = O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L
	integrand = lambda z : 1 / (1+z)     / sqrt(  O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L  )		

	#Integration
	z_obs = z
	z_cmb = 1089 #As Beta (not cmb). But 1089 (cmb) would be the exagerated maximun possible redshift for the birth 
	z_now = 0


	integral, error = quad( integrand , z_obs, z_cmb) #
	
	#t = ageoftheuniverse - (integral * (1 / H_sec) / secondsinyear)
	t = (integral * (1 / H_sec)) / secondsinyear

	return t



"""===================================================
Reddening functions        
==================================================="""


def BBB_nf2(bbb_x, bbb_y, BBebv, z ):

	"""
	
	## input:

	## output:

	"""
	#Application of reddening - reading E(B-V) from MCMC sampler
	RV= 2.72

	#converting freq to wavelenght, to be able to use prevots function instead on simple linera interpolation 
	redd_x =  2.998 * 1e10 / (10**(bbb_x)* 1e-8)
	redd_x= redd_x[::-1]

	#	Define prevots function for the reddenin law redd_k	
	def function_prevot(x, RV):
   		y=1.39*pow((pow(10.,-4.)*x),-1.2)-0.38 ;
   		return y 

	bbb_k = function_prevot(redd_x, RV)

	bbb_k= bbb_k[::-1]

	bbb_Lnu_red = bbb_y * 10**(-0.4 * bbb_k * BBebv)

	return bbb_x, bbb_Lnu_red
	
def GALAXY_nf2( gal_nu, gal_Fnu,GAebv):

	"""
	This function computes the effect of reddening in the galaxy template (Calzetti law)

	## input:
	-frequencies in log nu
	- Fluxes in Fnu
	- the reddening value E(B-V)_gal
	## output:

	"""

	RV = 4.05		
	wl = np.arange(0.122, 2.18, 0.04)
	redd_k=[]
	for i in range(len(wl)):
		if (wl[i]>0.12 and wl[i]<0.63):
			k =   2.659*(-2.156+(1.509/wl[i])-(0.198/(wl[i]**2))+(0.011/(wl[i]**3)) )+RV
		elif (wl[i]>0.63 and wl[i]<2.2):
			k =  2.659*(-1.857+(1.040/wl[i]))+RV
		redd_k.append(k)
	
	micron2cm = 1e-4
	redd_wl_rest = wl*micron2cm	
	redd_wl = redd_wl_rest 

	redd_f_r= 2.998 * 1e10 / (redd_wl)
	redd_f_r= 2.998 * 1e10 / (redd_wl)
	redd_f = np.log10(redd_f_r)[::-1]
	redd_k= np.array(redd_k)[::-1]


	reddening = interp1d(redd_f, redd_k, kind='linear', bounds_error=True)
	reddening2 = extrap1d(reddening)
	
	if (np.amax(gal_nu) - np.amax(redd_f)) <= 1e7:
		redd_x = gal_nu			
	else: 	
		redd_x = np.log10(gal_nu)	
	gal_k = reddening2(redd_x)

	gal_Fnu_red = gal_Fnu* 10**(-0.4 * gal_k * GAebv)
	

	return gal_nu, gal_Fnu_red




Angstrom = 1e10

def z2Dlum(z):

	"""
	Calculate luminosity distance from redshift.
	"""
	
	#Cosmo Constants
	
	O_m = 0.266
	O_r =  0.
	O_k= 0.
	O_L = 1. - O_m
	H_0 = 70. #km/s/Mpc
	H_sec = H_0 / 3.0857e19 
	# equation

	a = 1/(1+z)
	E = O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L
	integrand = lambda z : 1 / sqrt(O_m * (1+z)**3 + O_r *(1+z)**4 + O_k *(1+z) + O_L)	

	#integration

	z_obs = z
	z_now = 0

	c_cm = 2.997e10

	
	integral = quad( integrand , z_now, z_obs)	
	dlum_cm = (1+z)*c_cm/ H_sec * integral[0] 
	dlum_Mpc = dlum_cm/3.08567758e24

	return dlum_cm
   

"""---------------------------------------------
 			COMPUTED QUANTITIES
-----------------------------------------------"""

def stellar_info(chain, data):

	"""
	computes stellar masses and SFRs
	"""

	gal_do,  irlum_dict, nh_dict, BBebv_dict, SFRdict = data.dictkey_arrays #call dictionary info

	#relevanta parameters form the MCMC chain
	tau_mcmc = chain[:,0] 	
	age_mcmc = chain[:,1] 
	GA = chain[:,6] - 18. #1e18 is the common normalization factor used in parspace.ymodel in order to have comparable NORMfactors	

	z = data.z
	distance = z2Dlum(z)

	#constants
	solarlum = const.L_sun.to(u.erg/u.second) #3.839e33
	solarmass = const.M_sun

	Mstar_list=[]
	SFR_list=[]


	for i in range (len (tau_mcmc)):		
		N = 10**GA[i]* 4* pi* distance**2 / (solarlum.value)/ (1+z)

		gal_do.nearest_par2dict(tau_mcmc[i], 10**age_mcmc[i], 0.)
		tau_dct, age_dct, ebvg_dct=gal_do.t, gal_do.a,gal_do.e
		SFR_mcmc =SFRdict[tau_dct, age_dct]

		# Calculate Mstar. BC03 templates are normalized to M* = 1 M_sun. 
		# Thanks to Kenneth Duncan, and his python version of BC03, smpy
		Mstar = np.log10(N * 1) 
		#Calculate SFR. output is in [Msun/yr]. 
		SFR = N * SFR_mcmc
		SFR_list.append(SFR.value)	
		Mstar_list.append(Mstar)	

	return np.array(Mstar_list)	, np.array(SFR_list)


def stellar_info_array(chain_flat, data, Nthin_compute):

	"""
	computes arrays of stellar masses and SFRs
	"""

	Ns, Npar = np.shape(chain_flat)  
	chain_thinned = chain_flat[0:Ns:int(Ns/Nthin_compute),:]

	Mstar, SFR = stellar_info(chain_thinned, data)
	Mstar_list = []
	SFR_list = []

	for i in range(Nthin_compute):
		for j in range(int(Ns/Nthin_compute)):
			Mstar_list.append(Mstar[i])
			SFR_list.append(SFR[i])

	Mstar1 = np.array(Mstar_list)        
	SFR1 = np.array(SFR_list)
	return Mstar1, SFR1



def sfr_IR(logL_IR):
	#calculate SFR in solar M per year 

	#for an array ofluminosities
	if len(logL_IR)>1:
		SFR_IR_list =[]

		for i in range(len(logL_IR)):
			SFR = 3.88e-44* (10**logL_IR[i])
			SFR_IR_list.append(SFR)
		SFR_IR_array = np.array(SFR_IR_list)
		return SFR_IR_array
	#or just for one luminosity
	else:		
		SFR = 3.88e-44* (10**logL_IR)
		return SFR


"""---------------------------------------------
 	PROJECTION MODELS ON BAND FILTER CURVES
-----------------------------------------------"""


def filters1( model_nus, model_fluxes, filterdict, z ):	

	"""
	Projects the model SEDs into the filter curves of each photometric band.

	##input:
	- model_nus: template frequencies [log10(nu)]
	- model_fluxes: template fluxes [F_nu]
	- filterdict: dictionary with all band filter curves' information.
				  To change this, add one band and filter curve, etc,
				  look at DICTIONARIES_AGNfitter.py
	- z: redshift

	##output:
	- bands [log10(nu)]
	- Filtered fluxes at these bands [F_nu]
	"""

	bands, files_dict, lambdas_dict, factors_dict = filterdict
	filtered_model_Fnus = []


	# Costumize model frequencies and fluxes [F_nu]
	# to same units as filter curves (to wavelengths [angstrom] and F_lambda)
	model_lambdas = nu2lambda_angstrom(model_nus) * (1+z)
	model_lambdas =  model_lambdas[::-1]
	model_fluxes_nu =  model_fluxes[::-1]
	model_fluxes_lambda = fluxnu_2_fluxlambda(model_fluxes_nu, model_lambdas) 
	mod2filter_interpol = interp1d(model_lambdas, model_fluxes_lambda, kind = 'nearest', bounds_error=False, fill_value=0.)			

	# For filter curve at each band. 
	# (Vectorised integration was not possible -> different filter-curve-arrays' sizes)
 	for iband in bands:

		# Read filter curves info for each data point 
		# (wavelengths [angstrom] and factors [non])
		lambdas_filter = np.array(lambdas_dict[iband])
		factors_filter = np.array(factors_dict[iband])
		iband_angst = nu2lambda_angstrom(iband)

		# Interpolate the model fluxes to 
		#the exact wavelengths of filter curves
		modelfluxes_at_filterlambdas = mod2filter_interpol(lambdas_filter)
		# Compute the flux ratios, equivalent to the filtered fluxes: 
		# F = int(model)/int(filter)
		integral_model = trapz(modelfluxes_at_filterlambdas*factors_filter, x= lambdas_filter)
		integral_filter = trapz(factors_filter, x= lambdas_filter) 	
		filtered_modelF_lambda = (integral_model/integral_filter)

		# Convert all from lambda, F_lambda  to Fnu and nu	
		filtered_modelFnu_atfilter_i = fluxlambda_2_fluxnu(filtered_modelF_lambda, iband_angst)
		filtered_model_Fnus.append(filtered_modelFnu_atfilter_i)

	return bands, np.array(filtered_model_Fnus)


c=	2.997e8

def fluxlambda_2_fluxnu (flux_lambda, wl_angst):

	"""
	Calculate F_nu from F_lambda.
	"""
	flux_nu = flux_lambda * (wl_angst**2. ) / c /Angstrom
	return flux_nu


def fluxnu_2_fluxlambda (flux_nu, wl_angst):

	"""
	Calculate F_lambda from  F_nu.
	"""
	flux_lambda = flux_nu / wl_angst**2 *c * Angstrom

	return flux_lambda #in angstrom

def nu2lambda_angstrom(nus):

	"""
	Calculate wavelength [angstrom] from frequency [log Hz].
	"""

	lambdas = c / (10**nus) * Angstrom
	return lambdas

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:   
	    if ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])>0:
	        return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
            else:
		return 0	
	elif x > xs[-1]:
	    if ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])>0:
                return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
            else:
		return 0	
	else:
	    if interpolator(x)>0:
            	return interpolator(x)
	    else: return 0
    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike



