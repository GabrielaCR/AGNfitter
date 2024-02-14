AGNfitter 2.0 (AGNfitter-rx release)
========

**A Bayesian MCMC approach to fitting Spectral Energy Distributions of AGN and galaxies**

`(Calistro-Rivera et al. 2016, Martínez-Ramírez et al. 2024)`



**Welcome to AGNfitter!** 

AGNfitter is a Python algorithm implementing a Bayesian methods to fit the spectral energy distributions (SEDs) of active galactic nuclei (AGN) and galaxies from the **radio to X-rays**.
Through this method, you will be able to robustly disentangle the physical processes responsible for the emission of your sources.

You only need a catalog of photometric data (wavelengths, fluxes and errors), take a few decisions (if you wish), and you are ready to go (see Example). You can find detailed documentation soon here.

AGNfitter makes use of a large library of theoretical, empirical, and semi-empirical models to characterize both the nuclear and host galaxy emission simultaneously. The model consists of four physical emission components: the X-ray corona, an accretion disk, a torus of AGN heated dust, stellar populations, cold dust in star forming regions and synchroton emission from the AGN and star forming regions.  


Requirements
-------------

* Numpy 
* Matplotlib
* Scipy
* Astropy 
* emcee
* ultranest
* 

Installation
----------------

Installation can be done by cloning this Github repository.

After installation, let's do a quick test:

**1)** In `example/SETTINGS_AGNfitter.py`, go to `def CATALOG_settings()` and change 

    cat['path'] ='/Users/USER/AGNfitter/'
    
to your AGNfitter path. These test settings point to the example catalog contained in  `data/catalog_example.txt`.
    
**2)** In the terminal, go to your AGNfitter path  and start

    ./RUN_AGNfitter_multi.py  example/SETTINGS_AGNfitter.py
    
You should have a nice example in your `cat['path']/OUTPUT` folder. 

###
Either make sure that the root AGNfitter directory is on your `PATH` or specify the full path to `RUN_AGNfitter_multi.py`.
###

Quick start
------------

**TASK 0 (optional):** If you wish to have a working path other than the AGNfitter code path, please change 

    cat['workingpath'] = cat['path']
    
to your costumized working path.


**TASK 1:** In your working path, configure your settings creating a file `my_SETTINGS_AGNfitter.py`.
This file should be created based on the example in `example/SETTINGS_AGNfitter.py` (copy+paste).
To get AGNfitter running this is the ONLY file you need to modify.

*TASK 1a:* Specify your catalog's format in:

    def CATALOG_settings()
        cat['path'] ='/Users/USER/AGNfitter/'
        cat['filename'] = 'data/catalog_example.txt
        cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. 
        cat['name'] = 0#'ID'            ## If ASCII: Column index (int) of source IDs
        cat...

*TASK 1b:* To construct the dictionary  please go to

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
    
    
**TASK 2:** Run AGNfitter with

    ./RUN_AGNfitter_multi.py my_SETTINGS_AGNfitter.py
   
This will run AGNfitter in series. In general there are a few more runtime options (see below).

Done!  

Documentation
----------------
A careful documentation will be soon available. In the meantime, some notes are available in the [Wiki](https://github.com/GabrielaCR/AGNfitter/wiki) as response to questions asked by users.

Citation and License
----------------
Please cite `Calistro Rivera et al. (2016)` and `Martinez-Ramirez et al (2024, submitted)` if this code has achieved its purpose and contributed to your
research. 
The BibTeX entry for the paper is:

    @ARTICLE{2016ApJ...833...98C,
           author = {{Calistro Rivera}, Gabriela and {Lusso}, Elisabeta and {Hennawi}, Joseph F. and {Hogg}, David W.},
            title = "{AGNfitter: A Bayesian MCMC Approach to Fitting Spectral Energy Distributions of AGNs}",
          journal = {\apj},
         keywords = {galaxies: active, galaxies: nuclei, galaxies: statistics, methods: statistical, quasars: general, Astrophysics - Astrophysics of Galaxies, Astrophysics - Instrumentation and Methods for Astrophysics},
             year = 2016,
            month = dec,
           volume = {833},
           number = {1},
              eid = {98},
            pages = {98},
              doi = {10.3847/1538-4357/833/1/98},
    archivePrefix = {arXiv},
           eprint = {1606.05648},
     primaryClass = {astro-ph.GA},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2016ApJ...833...98C},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

AGNfitter is an open-source software made available under the MIT License. See
the LICENSE file for details.
