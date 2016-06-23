import numpy as np
import math

import time
import cPickle
from astropy import units as u 



"""
Funtions used by the administrator to create objects with model fluxes
from original files (Not included).
"""

def STARBURST_read (fn):

	#reading
	c = 2.997e10
	c_Angst = 3.34e-19 #(1/(c*Angstrom)
	
	dh_wl_rest, dh_Flambda =  np.loadtxt(fn, usecols=(0,1),unpack= True)
	dh_wl = dh_wl_rest 
	dh_nu_r = np.log10(c / (dh_wl * 1e-8)) 
	dh_Fnu = dh_Flambda * (dh_wl**2. )* c_Angst

	#reverse , in order to have increasing frequency
	dh_nus= dh_nu_r[::-1]
	dh_Fnu = dh_Fnu[::-1]

	return dh_nus, dh_Fnu

def TORUS_read(tor_file):

	distance= 1e27

	tor_nu_rest, tor_nuLnu = np.loadtxt(tor_file, skiprows=0, usecols=(0,1),unpack= True)
	tor_Lnu = tor_nuLnu / 10**(tor_nu_rest)	
	tor_Fnu = tor_Lnu /(4. * np.pi * distance**2.) 

	return tor_nu_rest, tor_Fnu 

def BBB_read(bbb_file):

	bbb_nu_log_rest, bbb_nuLnu_log = np.loadtxt(bbb_file, usecols=(0,1),unpack= True)
	bbb_nu_exp = 10**(bbb_nu_log_rest) 
	bbb_nu = np.log10(10**(bbb_nu_log_rest) )
	bbb_nuLnu= 10**(bbb_nuLnu_log)
	bbb_Lnu = bbb_nuLnu / bbb_nu

	bbb_x = bbb_nu
	bbb_y =	bbb_nuLnu  / bbb_nu_exp

	return bbb_x, bbb_y


def write_totalfiles_SB():
	path ='/Users/Gabriela/Desktop/AGNfitter_newcode/AGNfitter_oldversion/'
	path2 ='/Users/Gabriela/Desktop/AGNfitter/'

	irlum1, all_filenames1 = np.loadtxt(path+'models/STARBURST/DALE.list' , usecols=(2,0), dtype = 'S', unpack= True)
	irlum2, all_filenames2 = np.loadtxt(path+'models/STARBURST/CHARY_ELBAZ.list' , usecols=(2,0), dtype = 'S',unpack= True)


	irlum = np.concatenate((irlum1.astype(float), irlum2.astype(float)))
	all_filenames = np.concatenate((all_filenames1, all_filenames2))

	wave_array=[]
	sed_array=[]
	for i in range(len(all_filenames)):
		wave, SED = STARBURST_read(path+'models/STARBURST/'+all_filenames[i])
		wave_array.append(wave)
		sed_array.append(SED)

	wave = np.asarray(wave_array)
	SED = np.asarray(sed_array)
			

	np.savez(path2 + 'models/STARBURST/dalehelou_charyelbaz_files', irlum, wave, SED)
#	cPickle.dump(np.array([irlum, wave, SED]), f2, protocol=2)

	a= np.load(path2 + 'models/STARBURST/dalehelou_charyelbaz_files.npz')



def write_totalfiles_TOR():

	path ='/Users/Gabriela/Desktop/AGNfitter_newcode/AGNfitter_oldversion/'
	path2 ='/Users/Gabriela/Desktop/AGNfitter/'

	tor_list = path +'models/TORUS/torus_templates_list.dat'
	nh, all_filenames= np.loadtxt(tor_list , usecols=(0,1), unpack=True, dtype = ('S'))
	nh = nh.astype(float)

	wave_array=[]
	sed_array=[]
	for i in range(len(all_filenames)):
		wave, SED = TORUS_read(path+all_filenames[i])
		wave_array.append(wave)
		sed_array.append(SED)

	wave = np.asarray(wave_array)
	SED = np.asarray(sed_array)

	np.savez(path2 + 'models/TORUS/silva_v1_files', nh, wave, SED)


def write_totalfiles_BBB():
	path ='/Users/Gabriela/Desktop/AGNfitter_newcode/AGNfitter_oldversion/'
	path2 ='/Users/Gabriela/Desktop/AGNfitter/'

	filename= 	path + 'models/BBB/richardsbbb.dat'
	wave, SED = BBB_read(filename)
	

	np.savez(path2 + 'models/BBB/richards_files', wave, SED)




class MODEL:

	"""
	construct the objects for callin the model frequencies, fluxes nd parameters
	to construct the dictionaries.
	"""


	def __init__(self, physcomponent, path):
		self.component = physcomponent
		self.path =path

	def build(self):

		wave_array = []
		sed_array = []

		if self.component=='starburst':

			npzfile = np.load(self.path + 'models/STARBURST/dalehelou_charyelbaz_files.npz')
			irlum,sb_wave, sb_SED = npzfile['arr_0'], npzfile['arr_1'], npzfile['arr_2']

			self.wave = sb_wave

			self.SED = sb_SED
			self.irlum = irlum

		if self.component== 'galaxy':

			self.wave,self.tau, self.tg,self.SED, self.SFR = cPickle.load(file(self.path + 'models/GALAXY/bc03_275templates_files.pickle', 'rb')) 


		if self.component== 'torus':
			npzfile = np.load(self.path + 'models/TORUS/silva_v1_files.npz')

			nh, torus_wv, torus_SED = npzfile['arr_0'], npzfile['arr_1'], npzfile['arr_2']


			self.wave = torus_wv
			self.SED = torus_SED
			self.nh=nh

		if self.component== 'bbb':

			npzfile = np.load(self.path + 'models/BBB/richards_files.npz')
			bbb_wv, bbb_SED = npzfile['arr_0'], npzfile['arr_1']


			filename= 	self.path + 'models/BBB/richardsbbb.dat'
			self.wave=bbb_wv
			self.SED =bbb_SED




def construct(path):

	"""
	saves objects of class MODEL, into pickle files.
	"""

	SB = MODEL('starburst', path)
	SB.build()
	f = open(path + 'models/STARBURST/dalehelou_charyelbaz_v1.pickle', 'wb')
	cPickle.dump(SB, f, protocol=2)
	f.close()


	TO = MODEL('torus', path)
	TO.build()
	f1 = open(path + 'models/TORUS/silva_v1.pickle', 'wb')
	cPickle.dump(TO, f1, protocol=2)
	f1.close()

	BB = MODEL('bbb', path)
	BB.build()
	f2 = open(path + 'models/BBB/richards.pickle', 'wb')
	cPickle.dump(BB, f2, protocol=2)
	f2.close()


	GA = MODEL('galaxy', path)
	GA.build()
	f2 = open(path + 'models/GALAXY/bc03_275templates.pickle', 'wb')
	cPickle.dump(GA, f2, protocol=2)
	f2.close()




