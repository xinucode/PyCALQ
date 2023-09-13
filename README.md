# Meeting 09/13{#Meetings}
- Amy prefers that their is no .YML files, reading in hdf5 files is preffered. 
- Arjun: Understand how h5py files work.
- Look at Andre's files to get a general idea of how it works, and how to adapt to methods for Jupyter notebook.
- |-> Building in easy way to look at pivots and output pickle files, is this best format? Build something!
- What analysis is already there?
- Stability plots, pivot classes, and Data handling?
- Jo: start analyzing data for NN partial waves analysis, reproduce phase shift and spectrum plot.
- Andre has provided data, so lets recreate the plots.  
  -> S-wave analysis start to reproduce. Andre has stability plots produced. 

  
# Lüscher Scmüscher
Full analysis chain of the finite volume spectrum from two-point correlators to phase-shifts and other infinite-volume observables using the Lüscher formalism

Pipeline:
- first we pick operators to fit correlation functions
For NN analysis: There are 3 steps
  - single nucleon stability study
  - two-nucleon stability study
  - phase-shift analysis
  - 

Questions
- what is a data handler? I will look into that
- what does raw LQCD data look like, how is it handled?
- Is h5py as output to get LQCD data the best way? (down the road)

Ideas
