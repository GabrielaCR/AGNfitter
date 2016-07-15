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
* scipy
* Astropy version 1.2 or later (pip install --no-deps astropy)
* acor (pip install acor)

Installation
----------------

Installation can by cloning this Github repository

After installation, let's do a quick test:

**1)** Make a copy of `example/SETTINGS_AGNfitter.py`, go to `def CATALOG_settings()` and change 
    cat['path'] ='/Users/USER/Codes/AGNfitter/'
    cat['outpath'] = '/Users/USER/AGNfitter/' 
to your AGNfitter and chosen output paths respectively.
    
    
**2)** In terminal go to the AGNfitter folder and start

    RUN_AGNfitter_multi.py my_SETTINGS_AGNfitter.py
    
You should have a nice example in your `cat['outpath']/OUTPUT` folder. Either make sure that the root AGNfitter directory is on your `PATH` or specify the full path to `RUN_AGNfitter_multi.py`.

Quick start
------------

To get AGNfitter running the ONLY file you need to open and modify is **example/SETTINGS_AGNfitter.py**.
You need just one thing to start: the catalog of sources.

**TASK 1:** Configure your settings based on the example in `example/SETTINGS_AGNfitter.py`

*TASK 1a:* Specify your catalog's format in:

    def CATALOG_settings()
        cat['path'] ='/Users/USER/Codes/AGNfitter/'
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
    RUN_AGNfitter_multi.py my_SETTINGS_AGNfitter.py
this will run AGNfitter in series for the two sources in the example catalog. In general there are a few more runtime options. You can see them with `python RUN_AGNfitter_multi.py -h`:

              
                    XXXX
        ___________ XXX _________________________________________________
                    XX      
                    X     
                    X                       AGNfitter                     
                __ X __                    ---------                
            /**\   |   /**\                                          
        ... (*** =  o  = ***) ...........................................
            \**/__ | __\**/                                     
                    X              Fitting SEDs of AGN and Galaxies  
                    X             in a MCMC Approach 
                xx              (Calistro Rivera et al. 2016)    
                xx               
        _______ xxx______________________________________________________
            xxxx

        usage: RUN_AGNfitter_multi.py [-h] [-c NCPU] [-n SOURCENUMBER] [-i] [-o]
                                    AGNfitterSettings

        positional arguments:
        AGNfitterSettings     AGNfitter settings file

        optional arguments:
        -h, --help            show this help message and exit
        -c NCPU, --ncpu NCPU  number of cpus to use for multiprocessing
        -n SOURCENUMBER, --sourcenumber SOURCENUMBER
                                specify a single source number to run (this is the
                                line number in hte catalogue not the source id/name)
        -i, --independent     run independently per source, i.e. do not create a
                                global model dictionary
        -o, --overwrite       overwrite model files



To run AGNfitter for a *single source* in the catalog, specify the line number as the sourcenumber argument: e.g.
    RUN_AGNfitter_multi.py --sourcenumber 0 my_SETTINGS_AGNfitter.py

To run AGNfitter in *batch mode* using python's multiprocessing capability and improve the efficiency, run e.g on a machine with 8 cpu cores
    RUN_AGNfitter_multi.py --ncpu 8 my_SETTINGS_AGNfitter.py
    
To run AGNfitter in a *distributed mode* on a compute cluster with multiple machines and a queue system, e.g using a PBS array job to specify the calalog line numbers
    RUN_AGNfitter_multi.py --independent --sourcenumber $PBS_ARRAY_ID my_SETTINGS_AGNfitter.py
See the qsub example in `example/run_agnfitter.qsub`. Here the `--independent` flag is required so that each job produces it's own model dictionary at it's own redshift (i.e. each source does not recreated the model dictionaries for the entire catalog). This can be more efficient for large catalogs where the model dictionary creation (which is not paralellized) can take a long time.

Additionally, you can specify `--overwrite` if you wish to recreate any existing models dictionaries (in case you change the z arrays).

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
