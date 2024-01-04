import logging
import yaml
import os
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import tqdm


import sigmond
import fvspectrum.sigmond_util as sigmond_util
import fvspectrum.sigmond_operator_info.operator
import general.plotting_handler as ph
import fvspectrum.sigmond_data_handling.data_handler as data_handler
from fvspectrum.sigmond_data_handling.data_files import DataFiles 
from fvspectrum.sigmond_data_handling.correlator_data import CorrelatorData

doc = '''
average_corrs - a task a read in and automatically average over any Lattice QCD temporal correlator data files 
                given within the same irrep row and total momentum

input
--------------------
general:
  ensemble_id: cls21_c103               #required
  ensembles_file: /home/sarahski/latticeQCD/pycalq/pycalq.average_task/fvspectrum/sigmond_utils/ensembles.xml #generated by PyCALQ
  project_dir: /latticeQCD/raid3/sarahski/lqcd/C103_R005/test_pycalq_project #required
  sampling_info:                        #not required 
    mode: Jackknife                     #default Jackknife
  tweak_ensemble:                       #not required
    omissions: []                       #default []
    rebin: 1                            #default 1
average_corrs:                          #required
  raw_data_files:                       #required
  - /latticeQCD/raid3/ahanlon/data/cls21_c103/updated_stats/sigmond.fwd/cls21_c103/nucleon_S0.bin
  average_by_bins: false                #not required #default false
  average_hadron_irrep_info: true       #not required #default true
  average_hadron_spatial_info: true     #not required #default true
  bins_mode: true                       #not required #default true
  create_pdfs: true                     #not required #default true
  create_pickles: true                    #not required #default true
  create_summary: true                  #not required #default true
  erase_original_matrix_from_memory: false #not required #default false
  figheight: 6                          #not required #default 6
  figwidth: 8                           #not required #default 8
  generate_estimates: true              #not required #default true
  ignore_missing_correlators: true      #not required #default true
  plot: true                            #not required #default true
  tmax: 64                              #not required #default 64
  tmin: 0                               #not required #default 0
  write_data: true                      #not required #default true
'''

class SigmondAverageCorrs:
    @property
    def info(self):
        return doc
    
    @property
    def summary_file(self):
        return os.path.join(self.proj_dir_handler.plot_dir(), f"{self.task_name}_summary") #add channel? project name?
    
    def averaged_file( self, binned, channel = None ): #add average info
        if binned:
            subdir = 'bins'
        else:
            subdir = 'samples'
        if channel:
            return os.path.join(self.proj_dir_handler.data_dir(subdir),f"averaged_corrs.hdf5[{channel}]")
        else:
            return os.path.join(self.proj_dir_handler.data_dir(subdir),f"averaged_corrs.hdf5")
        
    def corr_data_file(self,corr):
        return os.path.join(self.proj_dir_handler.data_dir("estimates"), f"{corr}_correlator_estimates.csv")

    def effen_data_file(self,corr):
        return os.path.join(self.proj_dir_handler.data_dir("estimates"), f"{corr}_effenergy_estimates.csv")
    
    def corr_plot_file(self,corr, ptype):
        return os.path.join(self.proj_dir_handler.plot_dir(f"{ptype}s"), f"{corr}_correlator.{ptype}")

    def effen_plot_file(self,corr, ptype):
        return os.path.join(self.proj_dir_handler.plot_dir(f"{ptype}s"), f"{corr}_effenergy.{ptype}")

    def __init__( self, task_name, proj_dir_handler, general_configs, task_configs ):
        self.task_name = task_name
        self.proj_dir_handler= proj_dir_handler
        #initialize your task, store default input in self.proj_dir_handler.log_dir() (basically, throw the full possible input with all parameters where all the assumed parameters have been filled in in there)
        
        if not task_configs:
            logging.critical(f"No directory to view. Add 'raw_data_files' to '{task_name}' task parameters.")

        #check that raw_data_files are real files and not in project
        raw_data_files = []
        if 'raw_data_files' in task_configs.keys():
            raw_data_files = task_configs['raw_data_files']
        raw_data_files = sigmond_util.check_raw_data_files(raw_data_files, general_configs['project_dir'])

        #these wont change, all correlators can and should be considered hermetian, 
        #       and time separation is a cosmetic parameter
        self.hermitian = True
        self.time_separation = 1

        #these can and will matter but only for special cases. Will need extra care when coding up. 
        #only coding up if we come across an instance of needing such
        self.subtract_vev = False
        self.vev_const = 0.0
        self.effective_energy_type = 0 #0=TimeForward, 1=TimeSymmetric, 2=TimeBackward?

        self.project_info = sigmond_util.setup_project(general_configs,raw_data_files)
        #check that raw data files match what's in the data handler currently, get list of files in datahandler. 
        this_data_handler = data_handler.DataHandler(self.project_info)
        present_files = list(this_data_handler.raw_data_files.bl_corr_files)+list(this_data_handler.raw_data_files.bl_vev_files)
        present_files += list(this_data_handler.raw_data_files.bin_files)+list(this_data_handler.raw_data_files.sampling_files)
        if set(present_files)<=set(raw_data_files):
            needed_files = raw_data_files[:]
            [needed_files.remove(file) for file in present_files]

            data_files = DataFiles()
            for data_dir in needed_files:
                data_files += data_handler._find_data_files(data_dir)
            logging.info("Searching through found raw data files...")
            this_data_handler._raw_data += this_data_handler._findLaphData(data_files)
            if data_files.bin_files:
                logging.info("Reading bin data")
                for bin_file in tqdm.tqdm(data_files.bin_files):
                    this_data_handler._raw_data += this_data_handler._findSigmondData(bin_file, sigmond.FileType.Bins)

            if data_files.sampling_files:
                logging.info("Reading sampling data")
                for smp_file in tqdm.tqdm(data_files.sampling_files):
                    this_data_handler._raw_data += this_data_handler._findSigmondData(smp_file, sigmond.FileType.Samplings)

            this_data_handler.raw_data_files += data_files
        else:
            # del this_data_handler
            this_data_handler.project_info = self.project_info
            data_files = DataFiles()
            this_data_handler.raw_data_files = data_files
            this_data_handler._raw_data = CorrelatorData()
            this_data_handler._findRawData()

        
        #other params
        self.other_params = {
            'write_data': True,
            'create_pdfs': True,
            'create_pickles': True,
            'create_summary': True,
            'plot': True,
            'figwidth':8,
            'figheight':6,
            'average_hadron_spatial_info': True,
            'average_hadron_irrep_info': True,
            'bins_mode': True,
            'tmin':0,
            'tmax':64,
            'erase_original_matrix_from_memory': True,
            'ignore_missing_correlators': True,
            'generate_estimates': True,
            'average_by_bins': True,
        }
        sigmond_util.update_params(self.other_params,task_configs) #update other_params with task_params, 
                                                                        #otherwise fill in missing task params

        if not self.other_params['create_pdfs'] and not self.other_params['create_pickles'] and not self.other_params['create_summary']:
            self.other_params['plot'] = False
        
        #make yaml output
        logging.info(f"Full input written to '{os.path.join(proj_dir_handler.log_dir(), 'full_input.yml')}'.")
        with open( os.path.join(proj_dir_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
            yaml.dump({"general":general_configs, task_name: task_configs}, log_file)

    def run( self ):
        this_data_handler, mcobs_handler, mcobs_get_handler = sigmond_util.get_data_handlers(self.project_info)
        self.data_handler = this_data_handler
        averaged_channels = dict()
        log_output = dict()

        for channel in self.data_handler.raw_channels:
            if channel.is_averaged:
                logging.warning(f"Channel {str(channel)} is averaged already.")

            averaged_channel = channel.averaged
            if averaged_channel not in averaged_channels:
                averaged_channels[averaged_channel] = list()
                log_output[str(averaged_channel)] = list()
            
            if channel.is_averaged:
                log_output[str(averaged_channel)].append("WARNING: is averaged already.")
            averaged_channels[averaged_channel].append(channel)
            log_output[str(averaged_channel)].append(str(channel))

        log_path = os.path.join(self.proj_dir_handler.log_dir(), 'channels_combined_log.yml')
        logging.info(f"List of averaged channels written to '{log_path}'.")
        with open(log_path, 'w+') as log_file:
            yaml.dump(log_output, log_file)

        self.averaged_operators = {}
        log_output = dict()
        for avchannel, rawchannels in averaged_channels.items():
            if avchannel not in self.averaged_operators.keys():
                self.averaged_operators[avchannel] = {}
                log_output[str(avchannel)] = {}
            for rawchannel in rawchannels:
                operators = [op for op in self.data_handler.getChannelOperators(rawchannel)]
                ops_map = _getOperatorsMap(operators, avchannel, self.other_params['average_hadron_spatial_info'], 
                                           self.other_params['average_hadron_irrep_info'])
                for avop, rawops in ops_map.items():
                    if avop not in self.averaged_operators[avchannel].keys():
                        self.averaged_operators[avchannel][avop] = []
                        log_output[str(avchannel)][str(avop)] = list()
                    if type(rawops)==list:
                        self.averaged_operators[avchannel][avop] += rawops
                        log_output[str(avchannel)][str(avop)] += [str(op) for op in rawops]
                    else:
                        self.averaged_operators[avchannel][avop].append(rawops)
                        log_output[str(avchannel)][str(avop)].append(str(rawops))

        log_path = os.path.join(self.proj_dir_handler.log_dir(), 'operators_combined_log.yml')
        logging.info(f"List of averaged operators written to '{log_path}'.")
        with open(log_path, 'w+') as log_file:
            yaml.dump(log_output, log_file)

        save_to_self = not self.other_params['generate_estimates'] and self.other_params['plot']
        if save_to_self:
            self.data = {}

        file_created = False
        for avchannel in self.averaged_operators:
            if file_created:
                wmode = sigmond.WriteMode.Update
            else:
                wmode = sigmond.WriteMode.Overwrite
            result_ops = [op.operator_info for op in self.averaged_operators[avchannel]]
            input_ops = list()
            if save_to_self:
                self.data[avchannel] = {}
            for i in range(len(list(self.averaged_operators[avchannel].values())[0])):
                an_item = list()
                for op2 in self.averaged_operators[avchannel]:
                    an_item.append( (self.averaged_operators[avchannel][op2][i].operator_info,1.0))
                input_ops.append(an_item)

            if self.other_params['average_by_bins']:
                result_obs = sigmond.doCorrelatorMatrixSuperpositionByBins(mcobs_handler, input_ops, result_ops, self.hermitian, self.other_params['tmin'],
                                                            self.other_params['tmax'], self.other_params['erase_original_matrix_from_memory'],
                                                            self.other_params['ignore_missing_correlators'])
                decoy = sigmond.XMLHandler()
                mcobs_handler.writeBinsToFile(result_obs,self.averaged_file(self.other_params['average_by_bins'],repr(avchannel)),decoy,wmode, 'H')
            else:
                result_obs = sigmond.doCorrelatorMatrixSuperpositionBySamplings(mcobs_handler, input_ops, result_ops, self.hermitian, self.other_params['tmin'],
                                                            self.other_params['tmax'], self.other_params['erase_original_matrix_from_memory'],
                                                            self.other_params['ignore_missing_correlators'])
                decoy = sigmond.XMLHandler()
                mcobs_handler.writeSamplingValuesToFile(result_obs,self.averaged_file(self.other_params['average_by_bins'],repr(avchannel)),decoy,wmode, 'H')
            file_created = True
            
            if save_to_self or self.other_params['generate_estimates']:
                for avop1,avop2 in itertools.product(self.averaged_operators[avchannel],self.averaged_operators[avchannel]):
                    corr = sigmond.CorrelatorInfo(avop1.operator_info,avop2.operator_info)
                    corr_name = repr(corr).replace(" ","-")
                    estimates = sigmond.getCorrelatorEstimates(mcobs_handler,corr,self.hermitian,self.subtract_vev,sigmond.ComplexArg.RealPart, 
                                                                self.project_info.sampling_info.getSamplingMode())
                    if save_to_self:
                        self.data[channel][corr] = {}
                        self.data[channel][corr]["corr"] = sigmond_util.estimates_to_df(estimates)
                    else:
                        sigmond_util.estimates_to_csv(estimates, self.corr_data_file(corr_name) )
                    estimates = sigmond.getEffectiveEnergy(mcobs_handler,corr,self.hermitian,self.subtract_vev,sigmond.ComplexArg.RealPart, 
                                                            self.project_info.sampling_info.getSamplingMode(),self.time_separation,self.effective_energy_type,self.vev_const)
                    if save_to_self:
                        self.data[channel][corr]["effen"] = sigmond_util.estimates_to_df(estimates)
                    else:
                        sigmond_util.estimates_to_csv(estimates, self.effen_data_file(corr_name) )
        
        logging.info(f"Saved averaged correlators to {self.averaged_file(self.other_params['average_by_bins'])}.")
        if self.other_params['generate_estimates']:
            logging.info(f"Estimates generated in directory {self.proj_dir_handler.data_dir('estimates')}.")
            
    def plot( self ):
        #make plot for each correlator -> save to pickle and pdf
        if self.other_params['plot']:
            logging.info(f"Saving plots to directory {self.proj_dir_handler.plot_dir()}...")
        else:
            logging.info(f"No plots requested.")
            return
        
        plh = ph.PlottingHandler()
        if self.other_params['create_summary']:
            plh.create_summary_doc("Average Data")

        #set up fig object to reuse
        plh.create_fig(self.other_params['figwidth'], self.other_params['figheight'])
        
        #loop through same channels #make loading bar
        for channel in self.averaged_operators:
            if self.other_params['create_summary']:
                plh.append_section(str(channel))
            for avop1,avop2 in itertools.product(self.averaged_operators[channel],self.averaged_operators[channel]):
                corr = sigmond.CorrelatorInfo(avop1.operator_info,avop2.operator_info)
                corr_name = repr(corr).replace(" ","-")
                if not self.other_params['write_data'] and self.other_params['plot']:
                    df = self.data[channel][corr]["corr"]
                else:
                    df = pd.read_csv(self.corr_data_file(corr_name))

                plh.clf()
                plh.correlator_plot(df, 0, avop1, avop2) #0 for regular corr plot

                if self.other_params['create_pickles']:
                    plh.save_pickle(self.corr_plot_file( corr_name, "pickle"))
                if self.other_params['create_pdfs'] or self.other_params['create_summary']:
                    plh.save_pdf(self.corr_plot_file( corr_name, "pdf"))

                if not self.other_params['write_data'] and self.other_params['plot']:
                    df = self.data[channel][corr]["effen"]
                else:
                    df = pd.read_csv(self.effen_data_file(corr_name))

                plh.clf()
                plh.correlator_plot(df, 1, avop1, avop2) #1 for effective energy plot

                if self.other_params['create_pickles']:
                    plh.save_pickle(self.effen_plot_file( corr_name, "pickle"))
                if self.other_params['create_pdfs'] or self.other_params['create_summary']:
                    plh.save_pdf( self.effen_plot_file(corr_name, "pdf")) 

                if self.other_params['create_summary']:
                    plh.add_correlator_subsection(repr(corr),self.corr_plot_file( corr_name, "pdf"),
                                                    self.effen_plot_file( corr_name, "pdf"))

        if self.other_params['create_summary']:
            plh.compile_pdf(self.summary_file) 
         



#drew, stolen from drew. all below this, stolen from drew. 
def _getOperatorsMap(operators, averaged_channel, get_had_spat=False, get_had_irrep=False):
    op_map = dict()
    for operator in operators:
        averaged_op = _getAveragedOperator(operator, averaged_channel, get_had_spat, get_had_irrep)
        if averaged_op in op_map:
            logging.critical(f"Conflicting operators {operator} and {op_map[averaged_op]}")

        op_map[averaged_op] = operator

    return op_map

def _getAveragedOperator(operator, averaged_channel, get_had_spat=False, get_had_irrep=False):
    if operator.operator_type is sigmond.OpKind.GenIrrep:
        logging.critical("Averaging of GIOperators not currently supported")

    op_info = operator.operator_info.getBasicLapH()
    if op_info.getNumberOfHadrons() == 1:
        obs_name = f"{NAME_MAP[op_info.getFlavor()]}-{op_info.getHadronSpatialType(1)}_{op_info.getHadronSpatialIdNumber(1)}"
        obs_id = 0
    else:
        obs_name = ""
        for had_num in range(1, op_info.getNumberOfHadrons()+1):
            had_name = NAME_MAP[op_info.getHadronFlavor(had_num)]
            had_psq = op_info.getHadronXMomentum(had_num)**2 + op_info.getHadronYMomentum(had_num)**2 + op_info.getHadronZMomentum(had_num)**2
            had_str = str(had_psq)
            if get_had_spat:
                had_spat_type = op_info.getHadronSpatialType(had_num)
                had_spat_id = op_info.getHadronSpatialIdNumber(had_num)
                had_str += f"_{had_spat_type}{had_spat_id}"

            if get_had_irrep:
                had_irrep_str = op_info.getHadronLGIrrep(had_num)
                had_str += f"{had_irrep_str}" 

            obs_name += f"{had_name}({had_str})"

        obs_id = op_info.getLGClebschGordonIdNum()

    return fvspectrum.sigmond_operator_info.operator.Operator(averaged_channel.getGIOperator(obs_name, obs_id))


NAME_MAP = {
    'pion': 'P',
    'eta': 'e',
    'phi': 'p',
    'kaon': 'K',
    'kbar': 'k',
    'nucleon': 'N',
    'delta': 'D',
    'sigma': 'S',
    'lambda': 'L',
    'xi': 'X',
    'omega': 'O',
}