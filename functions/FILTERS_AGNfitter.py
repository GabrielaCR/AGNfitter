
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
import sys,os
import time
import cPickle
import shelve
from astropy import units as u 
from astropy.table import Table
from astropy.io import ascii



class FILTER:

	"""
	Class FILTER

	Used for creation of an object that contains the filter troughput and summarizes all information of a filter.
	Objects used only internally.

	##input: 
	- filtername 
	- filename
	- wavelength/frequency 
	- unit of wavelength/frequency 
	- description

	"""     

	def __init__(self, filtername, filename, freqwl_format, freqwl_unit, description):
		
		self.filtername = filtername
		self.ID = id(self)
		self.original_filename=filename
		self.description = description
		## Extract lambdas and throughput factors from files
		wl_or_freq, self.factors =  np.loadtxt(filename, usecols=(0,1),unpack= True)

		if freqwl_format == 'frequency':
		    newfilter_lambdas = wl_or_freq* freqwl_unit.to(u.Angstrom, equivalencies=u.spectral())
		elif freqwl_format == 'wavelength':
		    newfilter_lambdas = wl_or_freq* freqwl_unit.to(u.Angstrom)
		else:
		    print 'ERROR: Please enter the strings "frequency" or "wavelength" for the ADDfilters[freq/wl_format].'

		self.lambdas = newfilter_lambdas
		c= 2.997e8
		Angstrom = 1e10

		self.central_lambda = np.sum(self.lambdas*self.factors)/np.sum(self.factors)
		self.central_nu = float(np.log10((Angstrom*c)/self.central_lambda))


class FILTER_SET:
	"""
	Class FILTER_SET

	Object that contains a selected group of filter objects.

	##input: 
	- filterset name 
	- filter settings dictionary
	- the dictionary of all filter objects

	"""   
	def __init__(self, filterset_name, filtersdict, filters_objects_all):

		self.name = filterset_name
		if filtersdict['filterset'] == 'filterset_default':
			default_filters =['SPIRE500', 'SPIRE350', 'SPIRE250', 'MIPS24', 'IRAC4',\
		 'IRAC3', 'IRAC2', 'IRAC1', 'Ks_2mass', 'H_2mass', 'J_2mass',\
		 'Y_VISTA', 'z_SUBARU', 'i_SUBARU', 'r_SUBARU', 'B_SUBARU', 'u_CHFT',\
		  'GALEX_2500']

		 	filters_objects_chosen = [o for o in filters_objects_all.values() \
								  if o.filtername in default_filters]
		else:
			filters_objects_chosen = [o for o in filters_objects_all.values() \
								  if filtersdict[o.filtername]==True]
		self.filternames =[i.filtername for i in filters_objects_chosen]
		#dictionaries lambdas_dict, factors_dict
		lambdas_dict = defaultdict(list)
		factors_dict = defaultdict(list)
		central_nu_list=[]

		for i in range(len(filters_objects_chosen)):
			lambdas_dict[filters_objects_chosen[i].central_nu].append(filters_objects_chosen[i].lambdas)
			factors_dict[filters_objects_chosen[i].central_nu].append(filters_objects_chosen[i].factors)
			central_nu_list.append(filters_objects_chosen[i].central_nu)

		self.central_nu_array=np.array(sorted(central_nu_list))
		self.lambdas_dict= lambdas_dict
		self.factors_dict= factors_dict



	def save(self,filename):
		f = open(filename, 'wb')
		cPickle.dump(self, f, protocol=2)
		f.close()



def add_newfilters(filters_objects_all_filename, ADDfilters_dict, path): 
	"""
	Fuction called when new filters need to be added.

	##input:
	- (str) name of the dictionary containing all filter objects
	- dictionary from the settings file, giving all filters to add.
	- AGNfitter path.

	"""    
	## Create the file with ALL filters (no needed for the user as the file ALL_FILTERS is already provided )
	if not os.path.lexists(filters_objects_all_filename): 

		filters_objects_all=dict()

		for i in range(len(ADDfilters_dict['names'])):
			filtername=ADDfilters_dict['names'][i]
			filters_objects_all[filtername] = FILTER(filtername, ADDfilters_dict['filenames'][i], ADDfilters_dict['freq/wl_format'][i],\
												 ADDfilters_dict['freq/wl_unit'][i],ADDfilters_dict['description'][i])

		central_lambdas = [filters_objects_all[i].central_lambda for i in ADDfilters_dict['names']]
		central_nus = [filters_objects_all[i].central_nu for i in ADDfilters_dict['names']]
		IDs= [filters_objects_all[i].ID for i in ADDfilters_dict['names']]
		f = open(filters_objects_all_filename, 'wb')
		cPickle.dump(filters_objects_all, f, protocol=2) 

	## Add the user's new filters 
	else: 
		with open(filters_objects_all_filename, "rb") as f: ## Get old set of all filters
			filters_objects_all = cPickle.load(f)
		f.close()

		for i in range(len(ADDfilters_dict['names'])): ## Add the new ones
			if ADDfilters_dict['names'][i] not in filters_objects_all.keys():
				filtername=ADDfilters_dict['names'][i]
				filters_objects_all[filtername] = FILTER(filtername, ADDfilters_dict['filenames'][i], ADDfilters_dict['freq/wl_format'][i],\
													 ADDfilters_dict['freq/wl_unit'][i],ADDfilters_dict['description'][i])
			else:
				'ERROR in adding new filter: "'+ ADDfilters_dict['names'][i] + '" is already recorded.'	

	## Save the info of all filters, including new added ones in a table ALL_FILTERS_info.dat.
	filters_table = Table(map(list,zip(*[[o.ID, o.filtername, '{:.0f}'.format(o.central_lambda), '{:.3f}'.format(o.central_nu), o.description] for o in filters_objects_all.values()])) ,\
					names=('ID (disk)','filtername', 'central lambda (Angstrom)', 'central nu (log Hz)', 'description of filter') )
	filters_table_sorted =filters_table.sort('central nu (log Hz)')
	ascii.write(filters_table, path + 'models/FILTERS/ALL_FILTERS_info.dat', delimiter ='|', overwrite=True)
	
	## save again the dictionary of all FILTER objects, now including the new aded ones.
	a = open(filters_objects_all_filename, 'wb')
	cPickle.dump(filters_objects_all, a, protocol=2) ## save list of FILTER objects

	## Tell the user to add the new filters to the settings file,
	## So that these can be reused for other configurations.
	## And to rerun the code.
	filters_in_settings =["filters['"+i+"'] = True" for i in ADDfilters_dict['names']]
	filters_in_settings='\n'.join(filters_in_settings)
	sys.exit("FILTERS ADDED\n---------------\nYour filters "+ str(ADDfilters_dict['names'])+ "have been successfully added.\n>> Please COPY THE FOLOWING LINES in the list of filters in your setting file.\n\n"+ \
	  filters_in_settings+"\n\n" + \
	 ">> Change the setting to\n   filters_dict['add_filters']==False \n>> Run the code again. ")



def create_filtersets(filters_dict, path):
	"""
	Creates new filter_sets following the choice of the user given in the settings file.

	##input:
	- dictionary from the settings file, giving all filters to add.
	- AGNfitter path.

	## called in DICTIONARIES_AGNfitter.py

	"""   
	filters_objects_all_filename = path+filters_dict['path']+ 'ALL_FILTERS'
	filterspath = filters_dict['path']

	if not os.path.lexists(filters_objects_all_filename):
		ADDfilters_dict = filters_dict['add_filters_dict']
		add_newfilters(filters_objects_all_filename, ADDfilters_dict, path)

	if filters_dict['add_filters']==True:
		ADDfilters_dict = filters_dict['add_filters_dict']
		add_newfilters(filters_objects_all_filename, ADDfilters_dict, path)	    	

	filters_objects_all = cPickle.load(file(filters_objects_all_filename, 'rb'))
	filterset = FILTER_SET(filters_dict['filterset'], filters_dict, filters_objects_all)

	return filterset

		