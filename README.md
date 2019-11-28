AGNfitter 
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
