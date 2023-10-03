# pyCALC
Correlator Analysis and Luscher Quantization Condition.

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
usage: run.py [-h] [-g GENERAL] [-t TASKS [TASKS ...]]

options:
  -h, --help            show this help message and exit
  -g GENERAL, --general GENERAL
                        general configuration file
  -t TASKS [TASKS ...], --tasks TASKS [TASKS ...]
                        task(s) configuration file(s)
(base) PS C:\cygwin64\home\sarah\lqcd\luscher-scmuscher> run.py -g test
Hello World
general_configs: test
task_configs: None
(base) PS C:\cygwin64\home\sarah\lqcd\luscher-scmuscher> run.py -g test -t test
Hello World
general_configs: test
task_configs: ['test']
(base) PS C:\cygwin64\home\sarah\lqcd\luscher-scmuscher> run.py -g test -t test test2
Hello World
general_configs: test
task_configs: ['test', 'test2']
(base) PS C:\cygwin64\home\sarah\lqcd\luscher-scmuscher> ipython
Python 3.8.5 (default, Sep  3 2020, 21:29:08) [MSC v.1916 64 bit (AMD64)]
Type 'copyright', 'credits' or 'license' for more information
IPython 8.12.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: import luscher

In [2]: luscher.LuscherSchmusher( 'test' ).run()
---------------------------------------------------------------------------
AttributeError                            Traceback (most recent call last)
Cell In[2], line 1
----> 1 luscher.LuscherSchmusher( 'test' ).run()

AttributeError: module 'luscher' has no attribute 'LuscherSchmusher'

In [3]: luscher.LuscherSchmuscher( 'test' ).run()
Hello World
general_configs: test
task_configs: {}

In [4]: luscher.LuscherSchmuscher( 'test', 'test' ).run()
Hello World
general_configs: test
task_configs: test
```
