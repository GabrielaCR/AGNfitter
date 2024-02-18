Adding filters
=====

.. _installation:

Requirements 
------------



Adding photometric filters
----------------

.. note::
  Before adding the require filter please make sure that the filters you want to add are not already in the AGNfitter library. The filters are listed in the file >>file

To add new filters you can follow these steps:
* Find the transfer functions of the filters you want to add (e.g. here http://svo2.cab.inta-csic.es/theory/fps/ ) and check that the file has two columns: the first one for the wavelengths (must be in Angstroms) and the second one for the factors ( f (λ) with a maximum value of 1).
* Add the files of each filter in the main folder (/AGNfitter/).
* Make the following changes to the SETTINGS_AGNfitter.py file: filters[’add_filters’]= True, create a list with the
filters names in ADDfilters[’names’] y create a list with the filter file names in ADDfilters[’filenames’].
* Run the code with the following command: python RUN_AGNfitter_multi.py SETTINGS_AGNfitter.py. If it worked out well, you’ll get the following result:

>>> import 
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']


AAA

.. autofunction:: lumache.get_random_ingredients

AAA

.. autoexception:: lumache.InvalidKindError

For example:

>>> import 
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']
