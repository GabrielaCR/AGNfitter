AGNfitter
========
**A Bayesian MCMC approach to fitting Spectral Energy Distributions for AGN and galaxies**

Welcome to AGNfitter! 

AGNfitter is a Python algorithm implementing a fully Bayesian MCMC method to fit the spectral energy distributions (SEDs) of active galactic nuclei (AGN) and galaxies from the sub-mm to the UV.
Through this method, you will be able to robustly disentangle the physical processes responsible for the emission of your sources.

You only need a catalog of photometric data (wavelengths, fluxes and errors), take a few decisions (if you wish), and you are ready to go (see Example).

AGNfitter makes use of a large library of theoretical, empirical, and semi-empirical models to characterize both the nuclear and host galaxy emission simultaneously. The model consists of four physical emission components: an accretion disk, a torus of AGN heated dust, stellar populations, and cold dust in star forming regions. AGNfitter determines the posterior distributions of numerous parameters that govern the physics of AGN with a fully Bayesian treatment of errors and parameter degeneracies, allowing one to infer integrated luminosities, dust attenuation parameters, stellar masses, and star formation rates. A detailed presentation, testing and discussion on this can be found in `arxiv.link.`

Installation
----------------

Example
----------------
If you start with this code:

    @start AGNfitter
    
    

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
