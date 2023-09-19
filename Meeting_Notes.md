# Meeting 09/19 10 AM - Jo and Sarah
- Sarah had questions/comments about last weeks meeting notes
- 1. Sarah confirmed hdf5 is general data structure, but we need some type of input structure. Also make it simple enough that if you dont want to worry about setting up you dont need to.
- 2. For hdf5 files, important to use FAIR data (https://en.wikipedia.org/wiki/FAIR_data), METADATA derscriptors to note down important information for each fit. i.e. Original configuration file, irrep, quantum numbers, pivot info, rotation type, ... . Depends on the data handler/meta_keys being used. Relevant METADATA,  TODO 
- 2. Hdf5 files, look at constructing them and especially careful about closing them when using them, can lead to corruption. 
- Sarah is generally working on Sigmond_py_bindings. Emphasized these are not necessarily for our code, since we want to make it all on Python and fast, but might come in handy.
- Handler class, still a TODO. 
- How do we best save results fo reach file? Include a directory where ensembles are kept, then save an analysis run in a parent directory. Use pickle files to save in this parent directory, including a common structure for all analysis. TODO

## Sarah's todo list
- Understand Pickle
- Make simple python classes that I’ll think be useful for everyone and I will definitly use:
  - Config handler
  - Data Get Handler
  - Project Directory Handler
- Investigate how to install sigmond with pip
- Begin setting up sigmond tasks with pybindings


# Meeting 09/13 - Amy, Arjun, Jo
- 1. Amy prefers that their is no .YML files, reading in hdf5 files is preffered. 
- 2. Arjun: Understand how h5py files work. (any update?)
- Look at Andre's files to get a general idea of how it works, and how to adapt to methods for Jupyter notebook.
- |-> Building in easy way to look at pivots and output pickle files, is this best format? Build something!
- What analysis is already there?
- Stability plots, pivot classes, and Data handling?
- 3. Jo: start analyzing data for NN partial waves analysis, reproduce phase shift and spectrum plot.
- Andre has provided data, so lets recreate the plots.  
  -> S-wave analysis start to reproduce. Andre has stability plots produced.

  
