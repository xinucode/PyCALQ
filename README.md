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
- what is a data handler? I will look into that (Jo)
- what does raw LQCD data look like, how is it handled? Comes in hdf5 files usually, so data handler has to know about keys ...
- Is h5py as output to get LQCD data the best way? (down the road)

Ideas


Sample usage of proposed run setup:

```
(base) PS C:\cygwin64\home\sarah\lqcd\luscher-scmuscher> run.py -h
usage: run.py [-h] general_configs task_configs [task_configs ...]

positional arguments:
  general_configs  general configuration file
  task_configs     task(s) configuration file(s)

options:
  -h, --help       show this help message and exit
(base) PS C:\cygwin64\home\sarah\lqcd\luscher-scmuscher> run.py test test
Hello World
general_configs: test
task_configs: ['test']
(base) PS C:\cygwin64\home\sarah\lqcd\luscher-scmuscher> ipython
Python 3.8.5 (default, Sep  3 2020, 21:29:08) [MSC v.1916 64 bit (AMD64)]
Type 'copyright', 'credits' or 'license' for more information
IPython 8.12.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: import run

In [2]: run.LuscherSchmuscher( "test", "test").run()
Hello World
general_configs: test
task_configs: test
```