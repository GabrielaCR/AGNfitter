AGNfitter (beta_version - July 6)
========
**A Bayesian MCMC approach to fitting Spectral Energy Distributions of AGN and galaxies**

Welcome to AGNfitter! 

AGNfitter is a Python algorithm implementing a fully Bayesian MCMC method to fit the spectral energy distributions (SEDs) of active galactic nuclei (AGN) and galaxies from the sub-mm to the UV.
Through this method, you will be able to robustly disentangle the physical processes responsible for the emission of your sources.

You only need a catalog of photometric data (wavelengths, fluxes and errors), take a few decisions (if you wish), and you are ready to go (see Example).

AGNfitter makes use of a large library of theoretical, empirical, and semi-empirical models to characterize both the nuclear and host galaxy emission simultaneously. The model consists of four physical emission components: an accretion disk, a torus of AGN heated dust, stellar populations, and cold dust in star forming regions.  A detailed presentation, test and discussion of AGNfitter can be found in `https://arxiv.org/abs/1606.05648#`

Requirements
-------------
* Numpy version 1.6 or later
* Matplotlib version 1.4.0 or later
* Scipy
* Astropy version 1.2 or later (pip install --no-deps astropy)
* acor (pip install acor)

Installation
----------------

Installation can be done through Github
* Clone (recommended) 
* Downloading the Github tar (please, stay updated about changes and corrections, since this is still a beta version)

After installation, let's do a quick test:

**1)** Open **RUN_AGNfitter_multi.py**, go to def CATALOG_settings() and 
    change cat['path'] ='/Users/USER/Codes/AGNfitter/' to your path.
    
    
**2)** In terminal go to the AGNfitter folder and start

    ipython RUN_AGNfitter_multi.py
    
You should have a nice example in your OUTPUT folder.

Quick start
------------

To get AGNfitter running the ONLY file you need to open and modify is  **RUN_AGNfitter_multi.py**.
You need just one thing to start: the catalog of sources.

**TASK1:** Specify your catalog's format in:

    def CATALOG_settings()
        cat['path'] ='/Users/USER/Codes/AGNfitter/'
        cat['filename'] = 'data/catalog_example.txt
        cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. 
        cat['name'] = 0#'ID'            ## If ASCII: Column index (int) of source IDs
                                        ## If FITS: not yet
        cat...
**TASK2:** To construct the dictionary  please go to

    def FILTERS_settings():
        filters['dict_zarray'] = np.arange(zmin, zmax, zinterval)

Here you can specify the redshift ranges or a redshift array you need for you catalog.
The DICT_default only includes z=[0.283, 1.58] for the test. 
Please, consider this process takes around 0.1 minute per redshift element.
This process might be lengthy but you only have to do it once.

You can use the default combination of photometric bands by leaving

	filters['Bandset'] = 'BANDSET_default'.

Otherwise, if you like, you can specify the photometric bands included in your catalog by setting 

    def FILTERS_settings():
        ...
        filters['Bandset'] = 'BANDSET_settings' 
        
        filters['SPIRE500']= False
        filters['SPIRE350']= True        

and assigning 'True' to the keys corresponding to the photometric bands in your catalog.
    
**TASK3:** Decide if you want to fit only one source

    RUN_AGNfitter_onesource(sourceline)
    #RUN_AGNfitter_multiprocessing(nr. of processors)

or a total catalog, using many processors in a multi-core computer:

    #RUN_AGNfitter_onesource(sourceline)
    RUN_AGNfitter_multiprocessing(nr. of processors)
    
**TASK4:** In terminal go to the AGNfitter folder and start

    ipython RUN_AGNfitter_multi.py

Done!

Documentation
----------------
A careful documentation will be soon available.

Citation and License
----------------
Please cite `Calistro Rivera et al. (2016)`_ if this code has achieved its purpose and contributed to your
research. 
The BibTeX entry for the paper is:

    @ARTICLE{2016arXiv160605648C,
    author = {{Calistro Rivera}, G. and {Lusso}, E. and {Hennawi}, J.~F. and {Hogg}, D.~W.},
    title = "{AGNfitter: A Bayesian MCMC approach to fitting spectral energy distributions of AGN}",
    journal = {ArXiv e-prints},
    archivePrefix = "arXiv",
    eprint = {1606.05648},
    keywords = {Astrophysics - Astrophysics of Galaxies, Astrophysics - Instrumentation and Methods for Astrophysics},
    year = 2016,
    month = jun,
    adsurl = {http://adsabs.harvard.edu/abs/2016arXiv160605648C},
    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }



AGNfitter is an open-source software made available under the MIT License. See
the LICENSE file for details.
