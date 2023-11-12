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
        """
        Initializes the project directory handler.

        Parameters:
        - root (str): The root directory for the project.
        """
        self.root = root
        pathlib.Path(root).mkdir(parents=True, exist_ok=True) 

    def set_task(self, index, task_name ):
        """
        Sets up a subdirectory in the project directory for the given task.

        Parameters:
        - index (int): The index of the task.
        - task_name (str): The name of the task.

        Returns:
        - str: The path to the created task subdirectory.
        """
        subdir = os.path.join( self.root, f"{index}{task_name}")
        pathlib.Path(subdir).mkdir(parents=True, exist_ok=True) 
        self.task_name = task_name
        self.index = index
        self.current_working_dir = subdir
        return subdir
        
    def log_dir(self, sub_dir = []):
        """
        Sets up and returns the directory for log files.

        Parameters:
        - sub_dir (str or list): Additional subdirectories within the log directory.

        Returns:
        - str: The path to the created log directory.
        """
        if type(sub_dir) ==list:
            return self._set_task_subdir(["logs"]+sub_dir)
        elif type(sub_dir) ==str:
            return self._set_task_subdir(["logs",sub_dir])
        else:
            logging.error(f"{type(sub_dir)} is not a valid type for a subdirectory name in the log_dir of '{self.current_working_dir}'.")

    def data_dir(self, sub_dir = []):
        """
        Sets up and returns the directory for data files.

        Parameters:
        - sub_dir (str or list): Additional subdirectories within the data directory.

        Returns:
        - str: The path to the created data directory.
        """
        if type(sub_dir) ==list:
            return self._set_task_subdir(["data"]+sub_dir)
        elif type(sub_dir) ==str:
            return self._set_task_subdir(["data",sub_dir])
        else:
            logging.error(f"{type(sub_dir)} is not a valid type for a subdirectory name in the data_dir of '{self.current_working_dir}'.")

    def plot_dir(self, sub_dir = []):
        """
        Sets up and returns the directory for plot files.

        Parameters:
        - sub_dir (str or list): Additional subdirectories within the plot directory.

        Returns:
        - str: The path to the created plot directory.
        """
        if type(sub_dir) ==list:
            return self._set_task_subdir(["plots"]+sub_dir)
        elif type(sub_dir) ==str:
            return self._set_task_subdir(["plots",sub_dir])
        else:
            logging.error(f"{type(sub_dir)} is not a valid type for a subdirectory name in the plot_dir of '{self.current_working_dir}'.")

    def _set_task_subdir(self, sub_dirs):
        """
        Helper method to set up and return subdirectories within the current working directory.

        Parameters:
        - sub_dirs (list): List of subdirectory names.

        Returns:
        - str: The path to the created subdirectory.
        """
        dirs = [self.current_working_dir]+sub_dirs
        subdir = os.path.join(*dirs)
        pathlib.Path(subdir).mkdir(parents=True, exist_ok=True) 
        return subdir
    
    