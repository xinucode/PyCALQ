# PyCALQ
Correlator Analysis and Luscher Quantization Condition.

Full analysis chain of the finite volume spectrum from two-point correlators to phase-shifts and other infinite-volume observables using the Lüscher formalism

## Prerequisites

[sigmond pybindings (pip branch)](https://github.com/andrewhanlon/sigmond/tree/pip)

[sigmond scripts pybindings (pip branch)](https://github.com/andrewhanlon/sigmond_scripts/tree/pip)

## Setup
```
cd PyCALQ/
pip install -r requirements.txt
```

## Sample usage

```
(base) PS C:\cygwin64\home\sarah\lqcd\luscher-scmuscher> run.py -h
usage: run.py [-h] [-g GENERAL] [-t TASKS [TASKS ...]]

options:
  -h, --help            show this help message and exit
  -g GENERAL, --general GENERAL
                        general configuration file
  -t TASKS [TASKS ...], --tasks TASKS [TASKS ...]
                        task(s) configuration file(s)
```

## Tasks

general setup => general configs

### View Correlators
info

### Average Correlators
info

### Rotate Correlators
info

### Fit Correlators
info

### Compare Spectrums
info

### Single Channel Fit
info


## Skeleton Task Class
```
import logging

doc = '''
essential documentation
'''

class MyTaskName:
    @property
    def info(self):
        return doc

    def __init__( self, task_name, proj_files_handler, general_configs, task_configs ):
        self.proj_files_handler= proj_files_handler
        #initialize your task, store default input in self.proj_files_handler.log_dir() (basically, throw the full possible input with all parameters where all the assumed parameters have been filled in in there)

    def run( self ):
        pass
        # do the task, produce the data, data goes in self.proj_files_handler.data_dir(), info/warning/errors about the process goes in self.proj_files_handler.log_dir() (if any)

    def plot( self ):
        pass
        # make the plots, store in self.proj_files_handler.plot_dir(), again, any log/error warnings go in self.proj_files_handler.log_dir() as well (if any)
```
