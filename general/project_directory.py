import logging
import os
import pathlib

class ProjectDirectoryHandler:

    #Directory Structure:
    # /root
    #   /0task_name
    #     /logs
    #       /sub_dir
    #       /...
    #     /data
    #     /plots
    #   /1task_name
    #     /...

    def __init__(self, root ):
        self.root = root
        pathlib.Path(root).mkdir(parents=True, exist_ok=True) 

    #sets up a subdirectory in the project directory for the given task
    #internalizes the task subdirectory as 'current_working_dir'
    def set_task(self, index, task_name ):
        subdir = os.path.join( self.root, f"{index}{task_name}")
        pathlib.Path(subdir).mkdir(parents=True, exist_ok=True) 
        self.task_name = task_name
        self.index = index
        self.current_working_dir = subdir
        return subdir
        
    #sets up returns the 'log_dir' and any subdirectories given in string form
    #input -> sub_dir(str or list)
    #output -> creates directory root/0task_name/logs/sub_dir if given string
    #                        or root/0task_name/logs/sub_dir[0]/sub_dir[1]/... if given list
    def log_dir(self, sub_dir = []):
        if type(sub_dir) ==list:
            return self._set_task_subdir(["logs"]+sub_dir)
        elif type(sub_dir) ==str:
            return self._set_task_subdir(["logs",sub_dir])
        else:
            logging.error(f"{type(sub_dir)} is not a valid type for a subdirectory name in the log_dir of '{self.current_working_dir}'.")

    #sets up returns the 'data_dir' and any subdirectories given in string form
    #input -> sub_dir(str or list)
    #output -> creates directory root/0task_name/data/sub_dir if given string
    #                        or root/0task_name/data/sub_dir[0]/sub_dir[1]/... if given list
    def data_dir(self, sub_dir = []):
        if type(sub_dir) ==list:
            return self._set_task_subdir(["data"]+sub_dir)
        elif type(sub_dir) ==str:
            return self._set_task_subdir(["data",sub_dir])
        else:
            logging.error(f"{type(sub_dir)} is not a valid type for a subdirectory name in the data_dir of '{self.current_working_dir}'.")

    #sets up returns the 'plot_dir' and any subdirectories given in string form
    #input -> sub_dir(str or list)
    #output -> creates directory root/0task_name/plots/sub_dir if given string
    #                        or root/0task_name/plots/sub_dir[0]/sub_dir[1]/... if given list
    def plot_dir(self, sub_dir = []):
        if type(sub_dir) ==list:
            return self._set_task_subdir(["plots"]+sub_dir)
        elif type(sub_dir) ==str:
            return self._set_task_subdir(["plots",sub_dir])
        else:
            logging.error(f"{type(sub_dir)} is not a valid type for a subdirectory name in the plot_dir of '{self.current_working_dir}'.")

    def _set_task_subdir(self, sub_dirs):
        dirs = [self.current_working_dir]+sub_dirs
        subdir = os.path.join(*dirs)
        pathlib.Path(subdir).mkdir(parents=True, exist_ok=True) 
        return subdir
    
    