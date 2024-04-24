import logging
import os,glob
import pathlib

import general.task_manager as tm

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

    def __init__(self, root, tasks=[] ):
        """
        Initializes the project directory handler.

        Parameters:
        - root (str): The root directory for the project.
        - tasks (list): List of tasks to set up proactively
        """
        pathlib.Path(root).mkdir(parents=True, exist_ok=True) 

        self.root = root
        self.all_tasks = {}
        for i,task in enumerate(tasks):
            if os.path.exists(os.path.join(root,f"{i}{task}")):
                self.all_tasks[task] = ProjectDirectoryHandler(root)
                self.all_tasks[task].set_task(i,task)
                
    def set_task(self, index, task_name, recursive = True):
        """
        Sets up a subdirectory in the project directory for the given task.

        Parameters:
        - index (int): The index of the task.
        - task_name (str): The name of the task.
        - recursive (bool): Sets up tasks for self.all_tasks

        Returns:
        - str: The path to the created task subdirectory.
        """
        if task_name in self.all_tasks:
            self.task_name = task_name
            self.index = index
            self.current_working_dir = self.all_tasks[task_name].current_working_dir
            return self.all_tasks[task_name].current_working_dir
        else:
            subdir = os.path.join( self.root, f"{index}{task_name}")
            pathlib.Path(subdir).mkdir(parents=True, exist_ok=True) 
            self.task_name = task_name
            self.index = index
            self.current_working_dir = subdir
            if recursive:
                self.all_tasks[task_name] = ProjectDirectoryHandler(self.root)
                self.all_tasks[task_name].set_task(index, task_name, False)
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

    #configures the directory of potential outputs, designating two types of data as 'bins' or 'samples'
    def data_subdir( self, binned ):
        if binned:
            subdir = 'bins'
        else:
            subdir = 'samples'
        return self.data_dir(subdir)
    
    #Makes a consistent unique key for files
    def filekey(self, mom=None, rebin=1, sampling_type=None, rotate_type=None, tN=None, t0=None, tD=None, file_tag=""):
        mom_key = ""
        if mom!=None:
            mom_key = f"-PSQ{mom}"
        diag_key = ""
        if rotate_type!=None and tN!=None and t0!=None and tD!=None:
            diag_key = f"-{rotate_type}-{tN}tN-{t0}t0-{tD}tD"
        if file_tag:
            file_tag = "_"+file_tag
        return f"{self.task_name}{file_tag}-Nbin{rebin}{mom_key}{diag_key}_{sampling_type}"
    
    #filename for general samplings file
    def samplings_file( self, binned, channel = None, mom=None, rebin=1, sampling_type=None, rotate_type=None, tN=None, t0=None, tD=None, file_tag=""): #add rotation info, and then average info
        # if key:
        #     key += "_"
        if binned:
            sampling_type = 'bins'
        else:
            sampling_type = sampling_type+'-samplings' #'B' or 'J'
        basename = self.filekey(mom, rebin, sampling_type, rotate_type, tN, t0, tD, file_tag)+".hdf5"
        if channel:
            return os.path.join(self.data_subdir(binned),basename+f"[{channel}]")
        else:
            return os.path.join(self.data_subdir(binned),basename)
        
    #filename for csv file meant for computed estimates
    def estimates_file( self, key=""):
        if key:
            key += "_"
        return os.path.join(self.data_dir("estimates"),f"{self.task_name}_{key}estimates.csv")

    #file for correlator estimates
    def corr_estimates_file(self,corr):
        return os.path.join(self.data_dir("estimates"), f"{corr}_correlator_estimates.csv")

    #file for effective energy of correlator estimates
    def effen_estimates_file(self,corr):
        return os.path.join(self.data_dir("estimates"), f"{corr}_effenergy_estimates.csv")
    
    #filename for corr plot file
    def corr_plot_file(self,corr, ptype):
        return os.path.join(self.plot_dir(f"{ptype}s"), f"{corr}_correlator.{ptype}")

    #filename for effective energy plot file
    def effen_plot_file(self,corr, ptype):
        return os.path.join(self.plot_dir(f"{ptype}s"), f"{corr}_effenergy.{ptype}")
    
    #filename for correlator tmin plot file. Plots spectrum fit results for varying tmins
    def corr_fit_series_plot_file(self,corr,energy_type, ptype, series_type='tmin'):
        return os.path.join(self.plot_dir(f"{ptype}s"), f"{corr}_{energy_type}_{series_type}.{ptype}")

    #filename for latex pdf where all plots and relevant information is gathered
    def summary_file(self, mom=None):
        if mom!=None:
            return os.path.join(self.plot_dir(), f"{self.task_name}-PSQ{mom}_summary") #add channel? project name?
        else:
            return os.path.join(self.plot_dir(), f"{self.task_name}_summary") #add channel? project name?
    
    #filename for plot of the interacting energy levels of all spectrum fits grouped by irrep
    def summary_plot_file(self, ptype, filetag = ""):
        if filetag:
            filetag="-"+filetag
        return os.path.join(self.plot_dir(f"{ptype}s"), f"{self.task_name}{filetag}_summary_plot.{ptype}")
    
    #filename for important info for gevp pivot
    def pivot_file(self, rotate_type, tN, t0, tD, rebin=1, sampling_type=None, channel=None ):
        if tm.Task.rotate_corrs.name==self.task_name:
            return self.samplings_file(False, channel, None, rebin, sampling_type, rotate_type, tN, t0, tD, "pivot_info")
        else:
            return self.all_tasks[tm.Task.rotate_corrs.name].samplings_file(False, channel, None, rebin, sampling_type, rotate_type, tN, t0, tD, "pivot_info")
        
    #filename for histogram plot of operator overlaps
    def operator_overlaps_plot(self, op, ptype):
        return os.path.join(self.plot_dir(f"{ptype}s"), f"{op}_operator_overlaps.{ptype}")
    
    #filename for samplings of operator overlaps
    def operator_overlaps_samplings(self, channel = None, rebin=1, sampling_type=None, rotate_type=None, tN=None, t0=None, tD=None):
        if tm.Task.fit_spectrum.name==self.task_name:
            return self.samplings_file(False, channel, None, rebin, sampling_type, rotate_type, tN, t0, tD, "operator_overlaps")
        else:
            return self.all_tasks[tm.Task.fit_spectrum.name].samplings_file(False, channel, None, rebin, sampling_type, rotate_type, tN, t0, tD, "operator_overlaps")
    
    #finds all averaged correlator data files based on rebin, sampling_type == 'B' or 'J', and selects momentums in only_mom if desired
    def get_averaged_data(self,binned, rebin, sampling_type=None, only_mom = []):
        data_files = []
        if tm.Task.average_corrs.name==self.task_name:
            if only_mom:
                for mom in only_mom:
                    data_files += list(glob.glob(self.samplings_file(binned, None, mom, rebin,sampling_type)))
            else:
                data_files += list(glob.glob(self.samplings_file(binned, None, '*', rebin,sampling_type)))
            data_files.append(self.samplings_file(binned, None, None, rebin,sampling_type))
        else:
            if only_mom:
                for mom in only_mom:
                    data_files += list(glob.glob(self.all_tasks[tm.Task.average_corrs.name].samplings_file(binned, None, mom, rebin,sampling_type)))
            else:
                data_files += list(glob.glob(self.all_tasks[tm.Task.average_corrs.name].samplings_file(binned, None, '*', rebin,sampling_type)))
            data_files.append(self.all_tasks[tm.Task.average_corrs.name].samplings_file(binned, None, None, rebin,sampling_type))
        return data_files
    
    #finds all rotated correlator data files based on rebin, rotate_type = 'SP' or 'RP', 
        #diagonalization parameters tN, t0, and tD, and sampling_type = 'B' or 'J'.
    def get_rotated_data(self,binned, rebin, rotate_type, tN, t0, tD, sampling_type=None): #, only_mom = []
        data_files = []
        if tm.Task.rotate_corrs.name==self.task_name:
            data_files += list(glob.glob(self.samplings_file(binned, None, None, rebin,sampling_type, rotate_type, tN, t0, tD)))
        else:
            print(self.all_tasks[tm.Task.rotate_corrs.name].samplings_file(binned, None, None, rebin,sampling_type,rotate_type, tN, t0, tD))
            data_files += list(glob.glob(self.all_tasks[tm.Task.rotate_corrs.name].samplings_file(binned, None, None, rebin,sampling_type,rotate_type, tN, t0, tD)))
        return data_files
        
    #removes all old summary files for a given task
    def remove_summary_files(self):
        for f in glob.glob(self.summary_file('*')+".*"):
            os.remove(f)
        for f in glob.glob(self.summary_file()+".*"):
            os.remove(f)

    #removes all present averaged data files
    def remove_averaged_data(self,binned, rebin, sampling_type=None):
        if tm.Task.average_corrs.name==self.task_name:
            for f in glob.glob(self.samplings_file(binned, None, '*', rebin,sampling_type)):
                os.remove(f)
            for f in glob.glob(self.samplings_file(binned, None, '', rebin,sampling_type)):
                os.remove(f)
        else:
            for f in glob.glob(self.all_tasks[tm.Task.average_corrs.name].samplings_file(binned, None, '*', rebin,sampling_type)):
                os.remove(f)
            for f in glob.glob(self.all_tasks[tm.Task.average_corrs.name].samplings_file(binned, None, '', rebin,sampling_type)):
                os.remove(f)