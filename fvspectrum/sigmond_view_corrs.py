import sigmond
import numpy as np
import logging
from typing import NamedTuple
import general.sigmond_data_handling.data_handler as data_handler
import general.sigmond_data_handling.data_files
import os

class ProjectInfo(NamedTuple):
  project_dir: str
  raw_data_dirs: list
  ensembles_file: str
  echo_xml: bool
  bins_info: sigmond.MCBinsInfo
  sampling_info: sigmond.MCSamplingInfo
  data_files: general.sigmond_data_handling.data_files.DataFiles
  precompute: bool
  latex_compiler: str

class SigmondViewCorrs:

    #initialize
    def __init__(self, task_name, general_params, task_params, log_dir):
        self.task_name = task_name
        if not task_params:
            logging.critical("No directory to view. Add 'raw_data_files' to 'view_data' task parameters.")

        #check that raw_data_files are real files
        raw_data_files = task_params.pop('raw_data_files',[])
        if not raw_data_files:
            logging.critical("No directory to view. Add 'raw_data_files' to 'view_data' task parameters.")
        if type(raw_data_files)!=list:
            if os.path.isdir(raw_data_files) or os.path.isfile(raw_data_files):
                raw_data_files = [raw_data_files]
            else:
                logging.critical("Parameter 'raw_data_files' must be a real directory.")
        else:
            filtered_raw_data_files = list(filter( lambda file: os.path.isdir(file) or os.path.isfile(file),raw_data_files))
            if filtered_raw_data_files!=raw_data_files:
                logging.critical("Item in 'raw_data_files' must be a real files.")
            raw_data_files = filtered_raw_data_files

        #check that raw_data_files are not in project
        parent_path = os.path.abspath(general_params['project_dir'])
        for file in raw_data_files:
            child_path = os.path.abspath(file)
            if parent_path == os.path.commonpath([parent_path, child_path]):
                logging.critical(f"Data directory '{child_path}' cannot be a subdirectory of project directory '{parent_path}'")


        print(raw_data_files)

        #generate ensembles file?
        # project_info = ProjectInfo(
        #     project_dir=general_params['project_dir'], raw_data_dirs=raw_data_files, ensembles_file=ensembles_file,
        #     echo_xml=False, bins_info=bins_info, sampling_info=sampling_info, data_files=data_files,
        #     precompute=precompute, latex_compiler=latex_compiler)

        # self.data_handler = data_handler.DataHandler(project_info)

    def run(self, data_dir, log_dir):
        pass

    def plot(self, plot_dir, log_dir):
        pass
