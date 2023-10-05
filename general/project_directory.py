import logging
import os
import pathlib

class ProjectDirectoryHandler:

    def __init__(self, root ):
        self.root = root
        pathlib.Path(root).mkdir(parents=True, exist_ok=True) 

    def task_subdir(self, index, task_name, stage):
        subdir = os.path.join( self.root, f"{index}{task_name}", stage)
        pathlib.Path(subdir).mkdir(parents=True, exist_ok=True) 
        return subdir