import logging
import os
import h5py
import tqdm

import sigmond
import general.task_manager as tm
import fvspectrum.sigmond_util as sigmond_util
from sigmond_scripts import data_handler
from sigmond_scripts.data_files import DataFiles
from sigmond_scripts.correlator_data import CorrelatorData

# states the data dependence of all sigmond tasks
dependencies = {
    tm.Task.preview_corrs: [],
    tm.Task.average_corrs: [],
    tm.Task.rotate_corrs: [tm.Task.average_corrs],
    tm.Task.fit_spectrum: [tm.Task.average_corrs, tm.Task.rotate_corrs],
    tm.Task.toy_corrs: [tm.Task.average_corrs]
}
raw_data_dependence = [tm.Task.preview_corrs, tm.Task.average_corrs]

#three stored correlator types:
    # raw correlators
    # averaged correlators
    # rotated correlators

#search the files for sigmond data
def parse_data_files(data_list):
    loglevel = logging.getLogger().getEffectiveLevel()
    logging.getLogger().setLevel(logging.WARNING)
    data_files = DataFiles()
    for item in data_list:
        data_files += data_handler._find_data_files(item)
    logging.getLogger().setLevel(loglevel)
    return list(data_files.bl_corr_files)+list(data_files.bl_vev_files)+list(data_files.bin_files)+list(data_files.sampling_files)

#handler that deals with general sigmond info and also manages the data
    #that's stored in memory
class SigmondProjectHandler:

    index = 0

    def __init__(self, general_configs, tasks_present):

        #set up sigmond necessities
        self.project_info = sigmond_util.setup_project(general_configs)
        self.ensemble_info = sigmond_util.get_ensemble_info(general_configs)
        shared = True
        if "qcdcomp" in os.uname().nodename: #cmu computing cluster
            shared = False
        self.data_handler = data_handler.DataHandler(self.project_info, shared)
        self.tasks_present = tasks_present
        self.contains_raw_data = False
        self.contains_averaged_data = False
        self.contains_rotated_data = False

        default_nodes = nodes=os.cpu_count()
        if shared:
            default_nodes /= 2
            default_nodes = int(default_nodes)
            default_nodes = min(default_nodes,8)
        self.nodes = general_configs.pop('nnodes', default_nodes)
        general_configs['nnodes'] = self.nodes

        #these wont change, all correlators can and should be considered hermetian, 
        #       and time separation is a cosmetic parameter
        self.hermitian = True
        self.time_separation = 1

        #these can and will matter but only for special cases. Will need extra care when coding up. 
        #only coding up if we come across an instance of needing such
        self.subtract_vev = False
        self.vev_const = 0.0
        self.effective_energy_type = 0 #0=TimeForward, 1=TimeSymmetric, 2=TimeBackward?

    #for a given set of data, manage whether to keep, delete, or search for data
    def add_data(self, in_data_files, current_data_files, clear_func, handler_data_dir, find_func, corr_data):
        #files already loaded into memory
        present_files = list(current_data_files.bl_corr_files)+list(current_data_files.bl_vev_files)
        present_files += list(current_data_files.bin_files)+list(current_data_files.sampling_files) 

        #files desired for next task
        parsed_raw_data_files = parse_data_files(in_data_files)
        common_files = set(present_files).intersection(parsed_raw_data_files)

        #if no common files, clear memory and load new files
        if not common_files:
            if present_files:
                clear_func()
            handler_data_dir.clear()
            handler_data_dir += in_data_files
            find_func()
        else:
            #delete files that aren't needed
            if common_files<set(present_files): 
                remove_files = list(set(present_files[:]))
                [remove_files.remove(file) for file in common_files]
                for file in remove_files:
                    if file in current_data_files.bl_corr_files:
                        current_data_files._bl_corr_files -= set(file)
                    elif file in current_data_files.sampling_files:
                        current_data_files._sampling_files -= set(file)
                    elif file in current_data_files.bin_files:
                        current_data_files._bin_files -= set(file)
                    elif file in current_data_files.bl_vev_files:
                        current_data_files._bl_vev_files -= set(file)
            #search for and add files that are needed and not present
            if common_files<set(parsed_raw_data_files):
                needed_files = list(set(parsed_raw_data_files[:]))
                [needed_files.remove(file) for file in common_files]

                data_files = DataFiles()
                for data_dir in needed_files:
                    data_files += data_handler._find_data_files(data_dir)
                logging.info("Searching through found data files...")
                new_corr_data = self.data_handler._findLaphData(data_files)
                if data_files.bin_files:
                    logging.info("Reading bin data")
                    for bin_file in tqdm.tqdm(data_files.bin_files):
                        new_corr_data += self.data_handler._findSigmondData(bin_file, sigmond.FileType.Bins)

                if data_files.sampling_files:
                    logging.info("Reading sampling data")
                    for smp_file in tqdm.tqdm(data_files.sampling_files):
                        new_corr_data += self.data_handler._findSigmondData(smp_file, sigmond.FileType.Samplings)
        
                for channel in new_corr_data._data_files.keys():
                    if channel in corr_data._data_files:
                        corr_data._data_files[channel]+=new_corr_data._data_files[channel]
                    else:
                        corr_data._data_files[channel]=new_corr_data._data_files[channel]
                for corr_info, corr_data_info in new_corr_data._correlators.items():
                    corr_data.addCorrelator(corr_info, **corr_data_info._asdict())

                current_data_files._sampling_files = current_data_files._sampling_files.union(data_files._sampling_files)
                current_data_files._bin_files = current_data_files._bin_files.union(data_files._bin_files)
                current_data_files._bl_corr_files.update(data_files._bl_corr_files)
                current_data_files._bl_vev_files.update(data_files._bl_vev_files)

    #manages memory associated with 'raw' data and add whatever is in raw_data_files (list)
    def add_raw_data(self, raw_data_files):
        self.add_data(raw_data_files, self.data_handler.raw_data_files, self.clear_raw_data, self.data_handler.rel_raw_datadir,
                      self.data_handler._findRawData, self.data_handler._raw_data)
        self.contains_raw_data = True

    #removes datafiles associated with each channel within channel_list (list) from memory
    def remove_raw_data_channels(self, channel_list):
        for channel in channel_list:
            self.data_handler.removeRawDataChannel(channel)

    #clears all present raw data from memory
    def clear_raw_data(self):
        data_files = DataFiles()
        self.data_handler.raw_data_files = data_files
        self.data_handler._raw_data = CorrelatorData()
        self.contains_raw_data = False

    def add_averaged_data(self, data_files):
        self.add_data(data_files, self.data_handler.averaged_data_files, self.clear_averaged_data, self.data_handler.rel_averaged_datadir,
                      self.data_handler.findAveragedData, self.data_handler._averaged_data) 
        self.contains_averaged_data = True

    def clear_averaged_data(self):
        data_files = DataFiles()
        self.data_handler.averaged_data_files = data_files
        self.data_handler._averaged_data = CorrelatorData()
        self.contains_averaged_data = False

    def remove_averaged_data_channels(self, channel_list):
        for channel in channel_list:
            self.data_handler.removeAveragedDataChannel(channel)

    def add_rotated_data(self, data_files):
        self.add_data(data_files, self.data_handler.rotated_data_files, self.clear_rotated_data, self.data_handler.rel_rotated_datadir,
                      self.data_handler.findRotatedData,self.data_handler._rotated_data)
        self.contains_rotated_data = True

    def clear_rotated_data(self):
        data_files = DataFiles()
        self.data_handler.rotated_data_files = data_files
        self.data_handler._rotated_data = CorrelatorData()
        self.contains_rotated_data = False

    def remove_rotated_data_channels(self, channel_list):
        for channel in channel_list:
            self.data_handler.removeRotatedDataChannel(channel)

    #when switching tasks, check for dependency of future tasks, and clear data from memory otherwise
    def switch_tasks(self):
        if self.contains_raw_data:
            clear_raw = True
            for task in self.tasks_present[self.index+1:]:
                if task in raw_data_dependence:
                    clear_raw = False
            if clear_raw:
                self.clear_raw_data()
        if self.contains_averaged_data:
            if self.tasks_present[self.index] in dependencies:
                clear = True
                for task in self.tasks_present[self.index+1:]:
                    if tm.Task.average_corrs in dependencies[task]:
                        clear = False
                if clear:
                    self.clear_averaged_data()
            else:
                self.clear_averaged_data()
        if self.contains_rotated_data:
            if self.tasks_present[self.index] in dependencies:
                clear = True
                for task in self.tasks_present[self.index+1:]:
                    if tm.Task.rotate_corrs in dependencies[task]:
                        clear = False
                if clear:
                    self.clear_rotated_data()
            else:
                self.clear_rotated_data()
        self.index+=1
