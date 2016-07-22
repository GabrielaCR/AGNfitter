AGNfitter (beta_version - July 22)
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

    RUN_AGNfitter_multi.py  example/SETTINGS_AGNfitter.py
    
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

    RUN_AGNfitter_multi.py my_SETTINGS_AGNfitter.py
   
This will run AGNfitter in series. In general there are a few more runtime options (see below).

Done!  


Runtime options
------------

You can see them with `python RUN_AGNfitter_multi.py -h`:

              
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

**Single source**

To run AGNfitter for a *single source* in the catalog, specify the line number as the sourcenumber argument: e.g.

    RUN_AGNfitter_multi.py --sourcenumber 0 my_SETTINGS_AGNfitter.py

**Working on a computer with multiple cores**

To run AGNfitter in *batch mode* using python's multiprocessing capability and improve the efficiency, run e.g on a machine with 10 cpu cores

    RUN_AGNfitter_multi.py --ncpu 10 my_SETTINGS_AGNfitter.py
    
**Working on a computer cluster with multiple different machines**

To run AGNfitter in a *distributed mode* on a compute cluster with multiple machines, a shared filesystem and a queue system, e.g using a PBS array job to specify the calalog line numbers

    RUN_AGNfitter_multi.py --independent --sourcenumber $PBS_ARRAY_ID my_SETTINGS_AGNfitter.py
    
The `--independent` flag is required so that each job produces its own model dictionary at its own redshift (i.e. each machine does not recreate the model dictionaries for the entire catalog). This can be more efficient for large catalogs where the model dictionary creation (which is not paralellized) can take a long time.

Additionally, you can specify `--overwrite` if you wish to recreate any existing models dictionaries (in case you change the z arrays).

Model Construction options
------------

As seen above, one main difference among the runtime options is the construction of the dictionary of models at the different redshifts of the sources.
Since the dictionary construction may be a lengthy process for large catalogs (0.1 min per redshift), you can choose among some options.
To summarise, there are basically three ways this dictionary can be constructed:

**Dict with grid of redshifts**: In `filters['dict_zarray']` you can specify a grid of redshifts which roughly covers the redshift range of your catalog. This is recommended for very large catalogs or not accurate redshifts. Ideally, the grid cells should not be larger than the redshift uncertainty.

**Dict with array of redshifts**: In `filters['dict_zarray']` you can specify the exact array of the redshifts in your catalog. This is recommended for small catalogs or very accurate redshifts. 

**Dict independent** By choosing the option -i (--independent) single model dictionaries will be produced for each source independently at its own redshift.  These dictionaries will be stored in each source's folder. This is recommended for compute clusters with multiple machines. 

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
