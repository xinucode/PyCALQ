import logging
import os
import xmltodict
import xml.etree.ElementTree as ET
import matplotlib
import shutil
from typing import NamedTuple
import pandas as pd

import sigmond
import fvspectrum.sigmond_data_handling.data_files
import fvspectrum.sigmond_data_handling.data_handler as data_handler


class ProjectInfo(NamedTuple):
  project_dir: str
  raw_data_dirs: list
  ensembles_file: str
  echo_xml: bool
  bins_info: sigmond.MCBinsInfo
  sampling_info: sigmond.MCSamplingInfo
  data_files: fvspectrum.sigmond_data_handling.data_files.DataFiles
  precompute: bool
  latex_compiler: str

def check_raw_data_files(raw_data_files, project_dir):
    #check that raw_data_files are real files
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
    parent_path = os.path.abspath(project_dir)
    for file in raw_data_files:
        child_path = os.path.abspath(file)
        if parent_path == os.path.commonpath([parent_path, child_path]):
            logging.critical(f"Data directory '{child_path}' cannot be a subdirectory of project directory '{parent_path}'")

    return raw_data_files

def get_ensemble_info(general_params):
    #check for ensembles file
    ensemble_file_path = os.path.join( os.path.dirname(__file__), "sigmond_utils", "ensembles.xml" )
    if not os.path.isfile(ensemble_file_path):
        logging.error("Ensembles file cannot be found.")

    with open(ensemble_file_path, 'r') as f:
        ensembles = xmltodict.parse(f.read())
    ids = [item['Id'] for item in ensembles['KnownEnsembles']['Infos']['EnsembleInfo']]
    if general_params['ensemble_id'] not in ids:
        logging.critical(f"Ensemble Id not found, check your 'ensemble_id' parameter or add your ensemble info to '{ensemble_file_path}'.")
    general_params["ensembles_file"] = ensemble_file_path

    #get bins info
    ensemble_info = sigmond.MCEnsembleInfo(general_params['ensemble_id'], ensemble_file_path)
    return ensemble_info

def setup_project(general_params, raw_data_files = []):
    ensemble_info = get_ensemble_info(general_params)
    ensemble_file_path = general_params["ensembles_file"]
    
    if 'tweak_ensemble' in general_params.keys():
        bins_info_config = general_params['tweak_ensemble']
        if 'keep_first' in bins_info_config:
            new_bins_info = ET.Element('MCBinsInfo')
            new_bins_info.append(ensemble_info.xml())
            tweaks = ET.SubElement(new_bins_info,'TweakEnsemble')
            ET.SubElement(tweaks,'KeepFirst').text = str(bins_info_config['keep_first'])
            ET.SubElement(tweaks,'KeepLast').text = str(bins_info_config['keep_last'])
            bins_info = sigmond.MCBinsInfo(sigmond.XMLHandler().set_from_string(ET.tostring(new_bins_info)))
        else:
            bins_info = sigmond.MCBinsInfo(ensemble_info)
            bins_info.setRebin(bins_info_config.get('rebin', 1))
            bins_info.addOmissions(set(bins_info_config.get("omissions", [])))
    else:
        bins_info = sigmond.MCBinsInfo(ensemble_info)
        general_params['tweak_ensemble'] = {}
        general_params['tweak_ensemble']['rebin'] = 1
        general_params['tweak_ensemble']['omissions'] = []

    #get sampling info
    if 'sampling_info' in general_params.keys():
        sampling_info_config = general_params['sampling_info']
        try:
            sampling_mode = sigmond.SamplingMode.create(sampling_info_config['mode'].lower())
        except KeyError as err:
            logging.critical("Unknown sampling mode {}".format(err))


        logging.info(str(sampling_mode).replace(".", ": "))

    else:
        sampling_mode = sigmond.SamplingMode.Jackknife
        general_params['sampling_info'] = {}
        general_params['sampling_info']['mode'] = "Jackknife"

        
    if sampling_mode == sigmond.SamplingMode.Bootstrap:
        try:
            if 'seed' not in sampling_info_config.keys():
                sampling_info_config['seed'] = 0
            if 'boot_skip' not in sampling_info_config.keys():
                sampling_info_config['boot_skip'] = 0

            sampling_info = sigmond.MCSamplingInfo(
                sampling_info_config['number_resampling'], sampling_info_config['seed'],
                sampling_info_config['boot_skip'])

        except KeyError as err:
            logging.critical("Missing required key {}".format(err))
    else:
        sampling_info = sigmond.MCSamplingInfo()

    datafiles = fvspectrum.sigmond_data_handling.data_files.DataFiles()

    return ProjectInfo( project_dir=general_params['project_dir'], raw_data_dirs=raw_data_files, ensembles_file=ensemble_file_path,
            echo_xml=False, bins_info=bins_info, sampling_info=sampling_info, data_files=datafiles,
            precompute=True, latex_compiler=None)

def set_latex_in_plots(style_import):
    try:
        style_import.use(os.path.join( os.path.dirname(__file__), "spectrum_plotting_settings", "spectrum.mplstyle" ))
    except:
        logging.warning('Spectrum style file has been moved or corrupted. Please consider reinstalling PyCALQ.')
        
    if not shutil.which('latex'):
        matplotlib.rcParams['text.usetex'] = False
        logging.warning("Latex not found on system, please install latex to get nice looking matplotlib plots.")
        return False
    return True

def update_params( other_params, task_params):
    for param in other_params:
        if param in task_params:
            other_params[param] = task_params[param]
        else:
            task_params[param] = other_params[param]

def get_data_handlers(project_info):
    
    this_data_handler = data_handler.DataHandler(project_info)

    mcobs_handler_init = ET.Element("MCObservables")
    if this_data_handler.raw_data_files.bl_corr_files:
        bl_corr_xml = ET.SubElement(mcobs_handler_init,"BLCorrelatorData")
        #file list info
    if this_data_handler.raw_data_files.bl_vev_files:
        bl_vev_files_xml = ET.SubElement(mcobs_handler_init,"BLVEVData")
        #file list info
    if this_data_handler.raw_data_files.bin_files:
        bin_files_xml = ET.SubElement(mcobs_handler_init,"BinData")
        for filename in this_data_handler.raw_data_files.bin_files:
            ET.SubElement(bin_files_xml,"FileName").text = filename
    if this_data_handler.raw_data_files.sampling_files:
        sampling_files_xml = ET.SubElement(mcobs_handler_init,"SamplingData")
        for filename in this_data_handler.raw_data_files.sampling_files:
            ET.SubElement(sampling_files_xml,"FileName").text = filename

    mcobs_init = sigmond.XMLHandler()
    mcobs_init.set_from_string(ET.tostring(mcobs_handler_init).decode())
    mcobs_get_handler = sigmond.MCObsGetHandler(mcobs_init, project_info.bins_info,project_info.sampling_info)

    # filemap = sigmond.XMLHandler()
    # mcobs_get_handler.getFileMap(filemap)
    # print(filemap.output(0))

    mcobs_handler = sigmond.MCObsHandler(mcobs_get_handler, project_info.precompute)
    mcobs_handler.setSamplingBegin()

    return this_data_handler, mcobs_handler, mcobs_get_handler


def estimates_to_csv( estimates, data_file):
    df = estimates_to_df( estimates )
    df.to_csv(data_file, index=False, header=True)

def estimates_to_df( estimates ):
    pestimates = [ {
        "aTime": key,
        "FullEstimate": estimates[key].getFullEstimate(), #should MC estimate have upper error and lower error as well?
        "AverageEstimate": estimates[key].getAverageEstimate(),
        "SymmetricError": estimates[key].getSymmetricError(),
        "RelativeError": estimates[key].getRelativeError(),
        } for key in estimates ]
    df = pd.DataFrame.from_dict(pestimates) 
    return df