# import sigmond
import numpy as np
import logging
import general.sigmond_data_handling.data_handler as data_handler

class ProjectInfo(NamedTuple):
  project_dir: str
  raw_data_dirs: list
  ensembles_file: str
  echo_xml: bool
  bins_info: sigmond.MCBinsInfo
  sampling_info: sigmond.MCSamplingInfo
  data_files: data_handling.data_files.DataFiles
  precompute: bool
  latex_compiler: str

class SigmondViewCorrs:

    def __init__(self, task_name, general_params, task_params):
        self.task_name = task_name
        if not task_params:
            logging.critical("No directory to view. Add 'raw_data_dirs' to 'view_data' task parameters.")

        raw_data_dirs = task_params.pop('raw_data_dirs',None)
        if not raw_data_directory:
            logging.critical("No directory to view. Add 'raw_data_dirs' to 'view_data' task parameters.")

        #generate ensembles file?
        # project_info = ProjectInfo(
        #     project_dir=general_params['project_dir'], raw_data_dirs=raw_data_dirs, ensembles_file=ensembles_file,
        #     echo_xml=False, bins_info=bins_info, sampling_info=sampling_info, data_files=data_files,
        #     precompute=precompute, latex_compiler=latex_compiler)

        # self.data_handler = data_handler.DataHandler(project_info)

