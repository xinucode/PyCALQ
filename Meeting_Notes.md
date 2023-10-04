# Meeting 10/03 10 AM
- changed name to PyCALQ
- Sarah has began incoorporating some of Drew's code, but still a work in progress
- Jo is incoorporating Andre's code, but still a work in progress
- Sarah is trying to finish "preview_corr" task as a first task
- Jo is considering how to split his anaysis into 'tasks'
- Jo has begun setting up a data_get_handler

## Sarah to do
- review joseph's code and make comments maybe
- finish setting up sigmond tasks with pybindings -> view_data
- Make simple python classes that I’ll think be useful for everyone and I will definitly use:
  - Project Directory Handler
  - add checks on config inputs in ConfigHandler
- Investigate how to install sigmond with pip

# Meeting 09/26 10 AM
-  Set up Config handler for YML, XML, JSON.
-  Run general config file, and info for specific file.
-  ensemble, channel, nbins. Have ensemble. For any given channel, need generic information. Made tasks options with library argparse.
-  Logging library useful for tracking errors and debugging.  Logging.error and logging.warning
-  general_configs are information for ensemble. Tasks being checked by task_configs.
-  Drew has view data task, including the tasks.
-  Find the necessary tasks, and optional tasks. Defaults must be decided, SOON. General config is in config_handler.py
-  luscher_schmuscher.py, needed checks. Can process input and directs it to analyzing code.
-  Pickle files are useful for saving class objects to memory.
-  Comends on luscher

## Sarah's todo list
- review joseph's code and make comments maybe
- Begin setting up sigmond tasks with pybindings -> view_data
- Make simple python classes that I’ll think be useful for everyone and I will definitly use:
  - test config handler 
  - Data Get Handler
  - Project Directory Handler
- Investigate how to install sigmond with pip

# Meeting 09/19 10 AM - Jo and Sarah
- Sarah had questions/comments about last weeks meeting notes
- 1. Sarah confirmed hdf5 is general data structure, but we need some type of input structure. Also make it simple enough that if you dont want to worry about setting up you dont need to.
- 2. For hdf5 files, important to use FAIR data (https://en.wikipedia.org/wiki/FAIR_data), METADATA derscriptors to note down important information for each fit. i.e. Original configuration file, irrep, quantum numbers, pivot info, rotation type, ... . Depends on the data handler/meta_keys being used. Relevant METADATA,  TODO 
- 2. Hdf5 files, look at constructing them and especially careful about closing them when using them, can lead to corruption. 
- Sarah is generally working on Sigmond_py_bindings. Emphasized these are not necessarily for our code, since we want to make it all on Python and fast, but might come in handy.
- Handler class, still a TODO. 
- How do we best save results fo reach file? Include a directory where ensembles are kept, then save an analysis run in a parent directory. Use pickle files to save in this parent directory, including a common structure for all analysis. TODO

## Joseph's todo list
- Understand Andre's code to reproduce NN scattering analysis
- Upload Luscher code for single-channel scattering analysis, and change it for more general methods.
- Checkout Data HANDLER Operations
- FAIR data metakeys

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

  
