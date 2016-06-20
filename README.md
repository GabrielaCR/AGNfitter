AGNfitter
========
**A Bayesian MCMC approach to fitting Spectral Energy Distributions of AGN and galaxies**

Welcome to AGNfitter! 

AGNfitter is a Python algorithm implementing a fully Bayesian MCMC method to fit the spectral energy distributions (SEDs) of active galactic nuclei (AGN) and galaxies from the sub-mm to the UV.
Through this method, you will be able to robustly disentangle the physical processes responsible for the emission of your sources.

You only need a catalog of photometric data (wavelengths, fluxes and errors), take a few decisions (if you wish), and you are ready to go (see Example).

AGNfitter makes use of a large library of theoretical, empirical, and semi-empirical models to characterize both the nuclear and host galaxy emission simultaneously. The model consists of four physical emission components: an accretion disk, a torus of AGN heated dust, stellar populations, and cold dust in star forming regions.  A detailed presentation, test and discussion of AGNfitter can be found in `arxivlink`

Installation
----------------

Install just downloading the Github tar.

After installation, let's do a quick test:

**1)** go to def CATALOG_settings() and 
    change cat['path'] ='/Users/USER/Codes/AGNfitter/' to your path.
    
    
**2)** Do what is told in  **TASK3**,  of **Quick start**

You should have a nice example in your OUTPUT folder.

Quick start
------------

To get AGNfitter running the only file you need to open and modify is **RUN_AGNfitter_multi.py**.
You need just one thing to start: the catalog of sources.

**TASK1:** Specify your catalog's format in:

    def CATALOG_settings()
        cat['path'] ='/Users/USER/Codes/AGNfitter/'
        cat['filename'] = 'data/catalog_example.txt
        cat['filetype'] = 'ASCII' ## catalog file type: 'ASCII' or 'FITS'. 
        cat['name'] = 0#'ID'            ## If ASCII: Column index (int) of source IDs
                                        ## If FITS: not yet
        cat...
    
**TASK2:** Decide if you want to fit only one source

    RUN_AGNfitter_onesource(sourceline)
    #RUN_AGNfitter_multiprocessing(nr. of processors)

or a total catalog, using many processors:

    #RUN_AGNfitter_onesource(sourceline)
    RUN_AGNfitter_multiprocessing(nr. of processors)
    
**TASK3:** In terminal go to the AGNfitter folder and start

    ipython RUN_AGNfitter_multi.py

Done!

Documentation
----------------
A careful documentation will be soon available at my webpage.

Citation and License
----------------
Please cite `Calistro Rivera et al. (2016)`_ if this code has achieved its purpose and contributed to your
research. I would be also very thankful if you could enter your name and paper to our ` testimonials list`_.
The BibTeX entry for the paper is:

AGNfitter is an open-source software made available under the MIT License. See
the LICENSE file for details.
