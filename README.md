# PyCALQ
Correlator Analysis and Luscher Quantization Condition.

Full analysis chain of the finite volume spectrum from two-point correlators to phase-shifts and other infinite-volume observables using the Lüscher formalism

Pipeline:
- first we pick operators to fit correlation functions
For NN analysis: There are 3 steps
  - single nucleon stability study
  - two-nucleon stability study
  - phase-shift analysis

Questions
- what is a data handler? 
- Best tasks for both data handling and analysis
- Exact steps for analyzing correlator data?

## Prerequisites

[sigmond pybindings (pip branch)](https://github.com/andrewhanlon/sigmond/tree/pip)

## Setup
```
cd PyCALQ/
pip install -r requirements.txt
python setup.py
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
add more...
