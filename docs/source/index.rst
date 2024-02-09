Welcome to AGNfitter's documentation!
===================================

**AGNfitter** is a Python code for modelling the spectral energy distributions (SEDs) of active galactic nuclei (AGN) and galaxies.
While AGNfitter can be applied on the photometry of any galaxy, the code includes a state-of-the-art library of AGN physical models and is ideal to infer detailed physical properties for AGNs, such as quasars, Seyfert galaxies, etc. 

The AGNfitter-rx release
--------

The first version of AGNfitter (Calistro-Rivera et al. 2016) has been recently expanded to model radio-to-X-ray photometry in the AGNfitter-rx release (Martinez-Ramirez et al. 2024).

.. image:: SED_manyrealizations_Mrk876.jpg



Here is the SED-fitting output of the code for the radio-to-X-ray photometry of the nearby AGN Mrk876. The black markers represent the photometric data and the lines of different colours represent the best-fit model and uncertainties for the different physical components listed in the legend. Below we also show the residuals. Additionally the code provides the user with tables with the output values of the physical parameters for the AGN and host galaxy, as well as many other optional output information.




 Here a list of the main features of the code, highlighting the new addition in the AGNfitter-rx release.

- Radio-to-X-ray SED modelling *(new)*
- Easy inclusion of your favorite new models toallow model comparison *(new)*
- Easy inclusion of filters  *(new)*
- Easy inclusion of user-tailored priors  *(new)*
- Energy balance prior between the stellar population and IR dust emission components *(new)*
- Sampled Probability Density Functions (PDFs) for all physical parameters

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.

The AGNfitter-rx release
--------


Citing the code
--------
If AGNfitter is useful for your research project, please include a citation to Calistro-Rivera et al. (2016) *and* Martinez-Ramirez et al. (2024)  in any publications. 


Contents
--------
.. toctree::

   usage
   api
