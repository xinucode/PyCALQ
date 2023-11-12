import numpy as np
import logging
from typing import NamedTuple
import xmltodict
import os
import yaml
import xml.etree.ElementTree as ET
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import pathlib
import pickle
import shutil
import pylatex

import sigmond
import fvspectrum.sigmond_data_handling.data_handler as data_handler
import fvspectrum.sigmond_data_handling.data_files
import fvspectrum.spectrum_plotting_settings.settings as psettings
import fvspectrum.sigmond_utils.util as utils

doc = '''
preview_corrs - a task a read in and estimate/plot any Lattice QCD temporal correlator data files given

required inputs
-------------
general:
  ensemble_id: cls21_c103
  project_dir: /latticeQCD/raid3/sarahski/lqcd/C103_R005/test_pycalq_project
preview_corrs:
  raw_data_files:
  - /latticeQCD/raid3/ahanlon/data/cls21_c103/updated_stats/sigmond.fwd/cls21_c103/nucleon_S0.bin

optional inputs
-------------
general:
  ensemble_id: cls21_c103
  ensembles_file: /home/sarahski/latticeQCD/luscher-scmuscher/fvspectrum/sigmond_utils/ensembles.xml
  project_dir: /latticeQCD/raid3/sarahski/lqcd/C103_R005/test_pycalq_project
  sampling_info:
    mode: Jackknife
  tweak_ensemble:
    omissions: []
    rebin: 1
preview_corrs:
  create_pdfs: true
  create_pickles: true
  create_summary: true
  figheight: 6
  figwidth: 8
  info: true
  plot: true
  raw_data_files:
  - /latticeQCD/raid3/ahanlon/data/cls21_c103/updated_stats/sigmond.fwd/cls21_c103/nucleon_S0.bin
  write_data: true
'''

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
      
class SigmondPreviewCorrs:

    @property
    def info(self):
        return doc

    def corr_data_file(self,corr):
        return os.path.join(self.proj_handler.data_dir(), f"{corr}_correlator_estimates.csv")

    def effen_data_file(self,corr):
        return os.path.join(self.proj_handler.data_dir(), f"{corr}_effenergy_estimates.csv")
    
    def corr_plot_file(self,corr, ptype):
        return os.path.join(self.proj_handler.plot_dir(f"{ptype}s"), f"{corr}_correlator.{ptype}")

    def effen_plot_file(self,corr, ptype):
        return os.path.join(self.proj_handler.plot_dir(f"{ptype}s"), f"{corr}_effenergy.{ptype}")
    
    @property
    def summary_file(self):
        return os.path.join(self.proj_handler.plot_dir(), f"{self.task_name}_summary") #add channel? project name?
    
    #initialize
    def __init__(self, task_name, proj_handler, general_params, task_params):
        self.task_name = task_name
        self.proj_handler = proj_handler

        if not task_params:
            logging.critical(f"No directory to view. Add 'raw_data_files' to '{task_name}' task parameters.")

        #check that raw_data_files are real files and not in project
        raw_data_files = []
        if 'raw_data_files' in task_params.keys():
            raw_data_files = task_params['raw_data_files']
        raw_data_files = _check_raw_data_files(raw_data_files, general_params['project_dir'])

        #check for ensembles file
        ensemble_file_path = os.path.join( os.path.dirname(__file__), "sigmond_utils", "ensembles.xml" )
        if not os.path.isfile(ensemble_file_path):
            logging.error("Ensembles file cannot be found.")

        self.latex = True
        plt.style.use(os.path.join( os.path.dirname(__file__), "spectrum_plotting_settings", "spectrum.mplstyle" ))
        if not shutil.which('latex'):
            matplotlib.rcParams['text.usetex'] = False
            logging.warning("Latex not found on system, please install latex to get nice looking matplotlib plots.")
            self.latex = False

        #check that ensemble_id is in ensembles file and if not print out options or instruct user to add the relevant info to ensembles.xml
        with open(ensemble_file_path, 'r') as f:
            ensembles = xmltodict.parse(f.read())
        ids = [item['Id'] for item in ensembles['KnownEnsembles']['Infos']['EnsembleInfo']]
        if general_params['ensemble_id'] not in ids:
            logging.critical(f"Ensemble Id not found, check your 'ensemble_id' parameter or add your ensemble info to '{ensemble_file_path}'.")
        general_params["ensembles_file"] = ensemble_file_path

        #get bins info
        ensemble_info = sigmond.MCEnsembleInfo(general_params['ensemble_id'], ensemble_file_path)
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

        self.project_info = ProjectInfo(
            project_dir=general_params['project_dir'], raw_data_dirs=raw_data_files, ensembles_file=ensemble_file_path,
            echo_xml=False, bins_info=bins_info, sampling_info=sampling_info, data_files=datafiles,
            precompute=True, latex_compiler=None)
        
        #other params
        self.other_params = {
            'write_data': True,
            'create_pdfs': True,
            'create_pickles': True,
            'create_summary': True,
            'plot': True,
            'figwidth':8,
            'figheight':6,
        }
        for param in self.other_params:
            if param in task_params:
                self.other_params[param] = task_params[param]
            else:
                task_params[param] = self.other_params[param]

        if not self.other_params['create_pdfs'] and not self.other_params['create_pickles']:
            self.other_params['plot'] = False
        
        #make yaml output
        logging.info(f"Full input written to '{os.path.join(proj_handler.log_dir(), 'full_input.yml')}'.")
        with open( os.path.join(proj_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
            yaml.dump({"general":general_params, task_name: task_params}, log_file)

    def run(self):
        #actually read files now -> ssssslllloooowwwwww!!!!! maybe look into mpi/slurm jobs already
        self.data_handler = data_handler.DataHandler(self.project_info)

        mcobs_handler_init = ET.Element("MCObservables")
        if self.data_handler.raw_data_files.bl_corr_files:
            bl_corr_xml = ET.SubElement(mcobs_handler_init,"BLCorrelatorData")
            #file list info
        if self.data_handler.raw_data_files.bl_vev_files:
            bl_vev_files_xml = ET.SubElement(mcobs_handler_init,"BLVEVData")
            #file list info
        if self.data_handler.raw_data_files.bin_files:
            bin_files_xml = ET.SubElement(mcobs_handler_init,"BinData")
            for filename in self.data_handler.raw_data_files.bin_files:
                ET.SubElement(bin_files_xml,"FileName").text = filename
        if self.data_handler.raw_data_files.sampling_files:
            sampling_files_xml = ET.SubElement(mcobs_handler_init,"SamplingData")
            for filename in self.data_handler.raw_data_files.sampling_files:
                ET.SubElement(sampling_files_xml,"FileName").text = filename

        mcobs_init = sigmond.XMLHandler()
        mcobs_init.set_from_string(ET.tostring(mcobs_handler_init))
        mcobs_get_handler = sigmond.MCObsGetHandler(mcobs_init, self.project_info.bins_info,self.project_info.sampling_info)
        mcobs_handler = sigmond.MCObsHandler(mcobs_get_handler, self.project_info.precompute)

        ##get operators -> print to logfile
        log_path = os.path.join(self.proj_handler.log_dir(), 'ops_log.yml')
        ops_list = {}
        ops_list["channels"] = {str(channel):{"operators":[str(op) for op in self.data_handler.getChannelOperators(channel)]} for channel in self.data_handler.raw_channels }
        
        logging.info(f"Channels and operators list written to '{log_path}'.")
        with open(log_path, 'w+') as log_file:
            yaml.dump(ops_list, log_file)
        
        save_to_self = not self.other_params['write_data'] and self.other_params['plot']
        if save_to_self:
            self.data = {}
        
        if not self.other_params['write_data'] and not self.other_params['plot']:
            logging.warning("You have set 'write_data' to 'False' and 'plot' to 'False' thus making this task obsolete. Congrats.")
            return

        logging.info(f"Saving correlator estimates to directory {self.proj_handler.data_dir()}...")
        for channel in self.data_handler.raw_channels:
            if save_to_self:
                self.data[channel] = {}
            for op1 in self.data_handler.getChannelOperators(channel):
                if save_to_self:
                    self.data[channel][op1] = {}
                for op2 in self.data_handler.getChannelOperators(channel):
                    corr = sigmond.CorrelatorInfo(op1.operator_info,op2.operator_info)
                    corr_name = repr(corr).replace(" ","-")
                    #add hermitian and subvev as inputs but default as True and False
                    estimates = sigmond.getCorrelatorEstimates(mcobs_handler,corr,True,False,sigmond.ComplexArg.RealPart, self.project_info.sampling_info.getSamplingMode())
                    if save_to_self:
                        self.data[channel][op1][op2] = {}
                        self.data[channel][op1][op2]["corr"] = _estimates_to_df(estimates)
                    else:
                        _estimates_to_csv(estimates, self.corr_data_file(corr_name) )
                    #add timesep, eff energy type, and subtraction constant? as parameters
                    estimates = sigmond.getEffectiveEnergy(mcobs_handler,corr,True,False,sigmond.ComplexArg.RealPart, self.project_info.sampling_info.getSamplingMode(), 1,0,0.0)
                    if save_to_self:
                        self.data[channel][op1][op2]["effen"] = _estimates_to_df(estimates)
                    else:
                        _estimates_to_csv(estimates, self.effen_data_file(corr_name) )



    def plot(self):
        #make plot for each correlator -> save to pickle and pdf
        if self.other_params['plot']:
            logging.info(f"Saving plots to directory {self.proj_handler.plot_dir()}...")
        else:
            logging.info(f"No plots requested.")
            return
        
        if self.other_params['create_summary']:
            doc = utils.create_doc("Preview Data")

        #set up fig object to reuse
        fig = plt.figure(figsize=(self.other_params['figwidth'], self.other_params['figheight'])) #add inputs
        
        #loop through same channels #make loading bar
        for channel in self.data_handler.raw_channels:
            if self.other_params['create_summary']:
                doc.append(pylatex.Command("section",str(channel)))
            for op1 in self.data_handler.getChannelOperators(channel):
                for op2 in self.data_handler.getChannelOperators(channel):
                    corr = sigmond.CorrelatorInfo(op1.operator_info,op2.operator_info)
                    corr_name = repr(corr).replace(" ","-")
                    if not self.other_params['write_data'] and self.other_params['plot']:
                        df = self.data[channel][op1][op2]["corr"]
                    else:
                        df = pd.read_csv(self.corr_data_file(corr_name))

                    text_x = 0.3
                    text_y = 0.7
                    plt.clf()
                    plt.gcf().set_size_inches(self.other_params['figwidth'], self.other_params['figheight'])
                    plt.errorbar(x=df["aTime"],y=df["FullEstimate"],yerr=df["SymmetricError"], linewidth=0.0, elinewidth=1.5, capsize=5, color=psettings.colors[0], marker=psettings.markers[0] )
                    plt.ylabel(r"$C(t)$")
                    plt.xlabel(r"$t/a_t$")
                    plt.figtext(text_x,text_y+0.1,f"snk: {str(op1)}") #double check that I'm not messing up sink and source
                    plt.figtext(text_x,text_y,f"src: {str(op2)}")
                    plt.tight_layout()

                    if self.other_params['create_pickles']:
                        with open(self.corr_plot_file( corr_name, "pickle"), "wb") as f:
                            pickle.dump(fig, f)
                    #put each channel in a pdf
                    if self.other_params['create_pdfs'] or self.other_params['create_summary']:
                        plt.savefig( self.corr_plot_file( corr_name, "pdf"), transparent=True ) 

                    #add pickle and pdf as option -> if none are selected -> no plot

                    if not self.other_params['write_data'] and self.other_params['plot']:
                        df = self.data[channel][op1][op2]["effen"]
                    else:
                        df = pd.read_csv(self.effen_data_file(corr_name))

                    plt.clf()
                    plt.gcf().set_size_inches(self.other_params['figwidth'], self.other_params['figheight'])
                    plt.errorbar(x=df["aTime"],y=df["FullEstimate"],yerr=df["SymmetricError"], linewidth=0.0, elinewidth=1.5, capsize=5, color=psettings.colors[0], marker=psettings.markers[0] )
                    if self.latex:
                        plt.ylabel(r"$a_tE_{\textup{fit}}$")
                    else: 
                        plt.ylabel(r"$a_tE_{fit}$")
                    plt.xlabel(r"$t/a_t$")
                    plt.figtext(text_x,text_y+0.1,f"snk: {str(op1)}")
                    plt.figtext(text_x,text_y,f"src: {str(op2)}")
                    plt.tight_layout()

                    if self.other_params['create_pickles']:
                        with open(self.effen_plot_file( corr_name, "pickle"), "wb") as f:
                            pickle.dump(fig, f)
                    #put each channel in a pdf
                    if self.other_params['create_pdfs'] or self.other_params['create_summary']:
                        plt.savefig( self.effen_plot_file(corr_name, "pdf"), transparent=True ) 

                    if self.other_params['create_summary']:
                        with doc.create(pylatex.Subsection(repr(corr))):
                            with doc.create(pylatex.Figure(position='H')) as thisfig:
                                with doc.create(pylatex.SubFigure(
                                    position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as left_fig:
                                    left_fig.add_image(self.corr_plot_file( corr_name, "pdf"),width=pylatex.NoEscape(r'\linewidth'), placement=pylatex.NoEscape("\centering") )
                                with doc.create(pylatex.SubFigure(
                                    position='b', width=pylatex.NoEscape(r'0.5\linewidth'))) as right_fig:
                                    right_fig.add_image(self.effen_plot_file( corr_name, "pdf"),width=pylatex.NoEscape(r'\linewidth'), placement=pylatex.NoEscape("\centering") )

        if self.other_params['create_summary']:
            utils.compile_pdf(doc,self.summary_file) #detect compiler?
                        

def _check_raw_data_files(raw_data_files, project_dir):
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

def _estimates_to_csv( estimates, data_file):
    df = _estimates_to_df( estimates )
    df.to_csv(data_file, index=False, header=True)

def _estimates_to_df( estimates ):
    pestimates = [ {
        "aTime": key,
        "FullEstimate": estimates[key].getFullEstimate(), #should MC estimate have upper error and lower error as well?
        "AverageEstimate": estimates[key].getAverageEstimate(),
        "SymmetricError": estimates[key].getSymmetricError(),
        "RelativeError": estimates[key].getRelativeError(),
        } for key in estimates ]
    df = pd.DataFrame.from_dict(pestimates) 
    return df
