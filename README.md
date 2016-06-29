AGNfitter (beta_version)
========
**A Bayesian MCMC approach to fitting Spectral Energy Distributions of AGN and galaxies**

Welcome to AGNfitter! 

AGNfitter is a Python algorithm implementing a fully Bayesian MCMC method to fit the spectral energy distributions (SEDs) of active galactic nuclei (AGN) and galaxies from the sub-mm to the UV.
Through this method, you will be able to robustly disentangle the physical processes responsible for the emission of your sources.

You only need a catalog of photometric data (wavelengths, fluxes and errors), take a few decisions (if you wish), and you are ready to go (see Example).

AGNfitter makes use of a large library of theoretical, empirical, and semi-empirical models to characterize both the nuclear and host galaxy emission simultaneously. The model consists of four physical emission components: an accretion disk, a torus of AGN heated dust, stellar populations, and cold dust in star forming regions.  A detailed presentation, test and discussion of AGNfitter can be found in `https://arxiv.org/abs/1606.05648#`

Requirements
-------------
(A version with less dependencies will be published July 2016)

* Numpy version 1.6 or later
* Matplotlib version 1.4.0 or later
* scipy
* Astropy version 1.2 or later (pip install --no-deps astropy)
* acor (pip install acor)

Installation
----------------

Install just downloading the Github tar.

After installation, let's do a quick test:

**1)** Open **RUN_AGNfitter_multi.py**, go to def CATALOG_settings() and 
    change cat['path'] ='/Users/USER/Codes/AGNfitter/' to your path.
    
    
**2)** In terminal go to the AGNfitter folder and start

    ipython RUN_AGNfitter_multi.py
    
You should have a nice example in your OUTPUT folder.

Quick start
------------

To get AGNfitter running you need to open and modify the file **RUN_AGNfitter_multi.py**.
You need just one thing to start: the catalog of sources.

**TASK1:** Specify your catalog's format in:

    def CATALOG_settings()
        cat['path'] ='/Users/USER/Codes/AGNfitter/'
        cat['filename'] = 'data/catalog_example.txt
        cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. 
        cat['name'] = 0#'ID'            ## If ASCII: Column index (int) of source IDs
                                        ## If FITS: not yet
        cat...
**TASK2:** To construct the dictionary, for the redshift ranges you need for you data, please go to the file DICTIONARY_AGNfitter.py and to the CLASS MODELSDICT. Please change the line in the __Init__ function with

    self.z_array =np.arange(zmin, zmax, zinterval)

and specify the redshift ranges you need for you catalog. 
The DICT_default only includes z=[0.14,0.28] for the test. 
Please, consider this process takes around 1 minute per redshift element.
This processs might be lengthy but you only have to do it once.
    
**TASK3:** Decide if you want to fit only one source

    RUN_AGNfitter_onesource(sourceline)
    #RUN_AGNfitter_multiprocessing(nr. of processors)

or a total catalog, using many processors:

    #RUN_AGNfitter_onesource(sourceline)
    RUN_AGNfitter_multiprocessing(nr. of processors)
    
**TASK4:** In terminal go to the AGNfitter folder and start

    ipython RUN_AGNfitter_multi.py

Done!

Documentation
----------------
A careful documentation will be soon available at my webpage.

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
