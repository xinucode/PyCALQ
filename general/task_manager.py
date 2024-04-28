import os
from enum import Enum

class Task(Enum): #encode tasks into enum
    preview_corrs = 0
    average_corrs = 1
    rotate_corrs = 2
    fit_spectrum = 3
    toy_corrs = 4
    compare_spectrums = 5

