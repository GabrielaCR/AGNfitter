Welcome to AGNfitter's documentation!
===================================

**AGNfitter** is a Python code for modelling the spectral energy distributions (SEDs) of active galactic nuclei (AGN) and galaxies.
While AGNfitter can be applied on the photometry of any galaxy, the code includes a state-of-the-art library of AGN physical models and is ideal to infer detailed physical properties for quasars, Seyfert galaxies, blazars, and other AGNs. 
AGNfitter implements Bayesian methodologies, allowing to robustly recover the physical parameters driving the emission, as well as their uncertainties and degeneracies.
In this documentation you should find everything you need to get the code running for your sources.

The AGNfitter-rx release
--------
.. image:: SED_manyrealizations_Mrk876.jpg


The first version of AGNfitter (Calistro-Rivera et al 2016) has been recently expanded to model radio-to-X-ray photometry in the AGNfitter-rx release (Martinez-Ramirez et al. 2024). Here a list of the main features of the latest release of the code.

- Radio-to-X-ray modelling *(new)*
- Easy inclusion of your favorite new models toallow model comparison *(new)*
- Energy balance prior between the stellar population and IR dust emission components *(new)*
-
- 

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.

Citing the code
--------
If AGNfitter is useful for your research project, please include a citation to Calistro-Rivera et al. (2016) *and* Martinez-Ramirez et al. (2024)  in any publications. 


Contents
--------
.. toctree::

   usage
   api
