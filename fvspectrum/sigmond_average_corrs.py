import logging
import yaml
import os
import pandas as pd
import itertools
import tqdm
from pqdm.processes import pqdm


import sigmond
import general.plotting_handler as ph
import fvspectrum.sigmond_util as sigmond_util
from sigmond_scripts.analysis.operator_info import operator as operator_lib
from sigmond_scripts.analysis.data_handling import data_handler
from sigmond_scripts.analysis.data_handling.data_files import DataFiles
from sigmond_scripts.analysis.data_handling.correlator_data import CorrelatorData

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
'''

#to do: specify irreps to average
class SigmondAverageCorrs:
    @property
    def info(self):
        return doc
    
    #final samplings file for averaged data
    def averaged_file( self, binned, mom=None, channel = None ): 
        if self.project_handler.project_info.sampling_info.isJackknifeMode():
            sampling_mode = 'J'
        else:
            sampling_mode = 'B'
        return self.proj_file_handler.samplings_file(binned, channel, mom, 
                                                     self.project_handler.project_info.bins_info.getRebinFactor(),
                                                     sampling_mode)

    def __init__( self, task_name, proj_file_handler, general_configs, task_configs, sph ):
        #set up fundamental project objects
        self.task_name = task_name
        self.proj_file_handler= proj_file_handler
        self.project_handler = sph
        
        if not task_configs:
            logging.critical(f"No directory to view. Add 'raw_data_files' to '{task_name}' task parameters.")

        #check that raw_data_files are real files and not in project
        raw_data_files = []
        if 'raw_data_files' in task_configs.keys():
            raw_data_files = task_configs['raw_data_files']
        raw_data_files = sigmond_util.check_raw_data_files(raw_data_files, general_configs['project_dir'])

        #retrieve data
        self.project_handler.add_raw_data(raw_data_files)
        self.data_handler = self.project_handler.data_handler
        
        #other params
        self.other_params = {
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
            'task_tag': "",
            'separate_mom': False
        }
        sigmond_util.update_params(self.other_params,task_configs) #update other_params with task_params, 
                                                                        #otherwise fill in missing task params

        if not self.other_params['create_pdfs'] and not self.other_params['create_pickles'] and not self.other_params['create_summary']:
            self.other_params['plot'] = False

        #filter channels
        self.channels = self.data_handler.raw_channels[:]
        final_channels = sigmond_util.filter_channels(task_configs, self.channels)
        remove_channels = set(self.channels)-set(final_channels)
        self.project_handler.remove_raw_data_channels(remove_channels)
        self.channels = final_channels

        #make yaml output
        logging.info(f"Full input written to '{os.path.join(proj_file_handler.log_dir(), 'full_input.yml')}'.")
        with open( os.path.join(proj_file_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
            yaml.dump({"general":general_configs, task_name: task_configs}, log_file)


        # if self.other_params['task_tag']:
        #     self.other_params['task_tag'] = "-"+self.other_params['task_tag']

    def run( self ):
        #get sigmond data handler
        mcobs_handler, mcobs_get_handler = sigmond_util.get_mcobs_handlers(self.data_handler,self.project_handler.project_info)
        averaged_channels = dict()
        log_output = dict()

        #check for already averaged channels, split channels in sets to be averaged
        for channel in self.channels:
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

        #put channel averaged info into log 
        log_path = os.path.join(self.proj_file_handler.log_dir(), 'channels_combined_log.yml')
        logging.info(f"List of averaged channels written to '{log_path}'.")
        with open(log_path, 'w+') as log_file:
            yaml.dump(log_output, log_file)

        #generate the lists of operators to average together
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

        #put list of operators averaged into logfile
        log_path = os.path.join(self.proj_file_handler.log_dir(), 'operators_combined_log.yml')
        logging.info(f"List of averaged operators written to '{log_path}'.")
        with open(log_path, 'w+') as log_file:
            yaml.dump(log_output, log_file)

        #determine if save data internally
        save_to_self = not self.other_params['generate_estimates'] and self.other_params['plot']
        if save_to_self:
            self.data = {}

        #for organizing momentum separated data files
        if self.other_params['separate_mom']:
            file_created = [False, False, False, False, False, False, False, False, False, False]
        else:
            file_created = [False]


        #do the averaging
        index = 0
        mom_key = None
        self.moms = []
        if self.other_params['separate_mom']:
            logging.info(f"Saving averaged correlators to {self.averaged_file(self.other_params['average_by_bins'],'*')}...")
        else:
            logging.info(f"Saving averaged correlators to {self.averaged_file(self.other_params['average_by_bins'])}...")
        for avchannel in tqdm.tqdm(self.averaged_operators):
            if self.other_params['separate_mom']:
                index = avchannel.momentum_squared
                mom_key = index
                self.moms.append(index)
                # self.other_params['task_tag'] = self.other_params['task_tag'][:-1]+str(index)
            if file_created[index]:
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
                result_obs = sigmond.doCorrelatorMatrixSuperpositionByBins(mcobs_handler, input_ops, result_ops, self.project_handler.hermitian, self.other_params['tmin'],
                                                            self.other_params['tmax'], self.other_params['erase_original_matrix_from_memory'],
                                                            self.other_params['ignore_missing_correlators'])
                decoy = sigmond.XMLHandler()
                mcobs_handler.writeBinsToFile(result_obs,self.averaged_file(self.other_params['average_by_bins'],mom_key,repr(avchannel)),decoy,wmode, 'H')
            else:
                result_obs = sigmond.doCorrelatorMatrixSuperpositionBySamplings(mcobs_handler, input_ops, result_ops, self.project_handler.hermitian, self.other_params['tmin'],
                                                            self.other_params['tmax'], self.other_params['erase_original_matrix_from_memory'],
                                                            self.other_params['ignore_missing_correlators'])
                decoy = sigmond.XMLHandler()
                mcobs_handler.writeSamplingValuesToFile(result_obs,self.averaged_file(self.other_params['average_by_bins'],mom_key,repr(avchannel)),decoy,wmode, 'H')
            file_created[index] = True
            
            #generate estimates for writing to file or plotting
            if save_to_self or self.other_params['generate_estimates']:
                for avop1,avop2 in itertools.product(self.averaged_operators[avchannel],self.averaged_operators[avchannel]):
                    corr = sigmond.CorrelatorInfo(avop1.operator_info,avop2.operator_info)
                    corr_name = repr(corr).replace(" ","-")
                    estimates = sigmond.getCorrelatorEstimates(mcobs_handler,corr,self.project_handler.hermitian,self.project_handler.subtract_vev,
                                                               sigmond.ComplexArg.RealPart, 
                                                                self.project_handler.project_info.sampling_info.getSamplingMode())
                    if save_to_self:
                        if avop1 not in self.data[avchannel]:
                            self.data[avchannel][avop1] = {}
                        self.data[avchannel][avop1][avop2] = {}
                        self.data[avchannel][avop1][avop2]["corr"] = sigmond_util.estimates_to_df(estimates)
                    else:
                        sigmond_util.estimates_to_csv(estimates, self.proj_file_handler.corr_estimates_file(corr_name) )
                    estimates = sigmond.getEffectiveEnergy(mcobs_handler,corr,self.project_handler.hermitian,self.project_handler.subtract_vev,
                                                           sigmond.ComplexArg.RealPart, self.project_handler.project_info.sampling_info.getSamplingMode(),
                                                           self.project_handler.time_separation,self.project_handler.effective_energy_type,
                                                           self.project_handler.vev_const)
                    if save_to_self:
                        self.data[avchannel][avop1][avop2]["effen"] = sigmond_util.estimates_to_df(estimates)
                    else:
                        sigmond_util.estimates_to_csv(estimates, self.proj_file_handler.effen_estimates_file(corr_name) )
        
        logging.info(f"Saved averaged correlators to {self.averaged_file(self.other_params['average_by_bins'])}.")
        if self.other_params['generate_estimates']:
            logging.info(f"Estimates generated in directory {self.proj_file_handler.data_dir('estimates')}.")
            
    def plot( self ):
        #determine if plots requested
        if self.other_params['plot']:
            logging.info(f"Saving plots to directory {self.proj_file_handler.plot_dir()}...")
        else:
            logging.info(f"No plots requested.")
            return
        
        #set up plotting handler
        plh = ph.PlottingHandler()
        self.moms = list(set(self.moms))

        #set up fig object to reuse
        plh.create_fig(self.other_params['figwidth'], self.other_params['figheight'])
        
        #loop through same channels #make loading bar
        if not self.other_params['generate_estimates'] and self.other_params['plot']:
            for channel in tqdm.tqdm(self.averaged_operators):
                self.write_channel_plots_data(channel, plh)
        else:
            for channel in tqdm.tqdm(self.averaged_operators):
                self.write_channel_plots(channel, plh) 
        
        #put all plots into summary doc
        if self.other_params['create_summary']:
            if self.other_params['separate_mom']:
                for i in self.moms:
                    plh.create_summary_doc("Average Data")
            else:
                plh.create_summary_doc("Average Data")
        
            index = 0
            for channel in self.averaged_operators:
                index = self.moms.index(channel.momentum_squared)
                plh.append_section(str(channel), index)
                for avop1,avop2 in itertools.product(self.averaged_operators[channel],self.averaged_operators[channel]):
                    corr = sigmond.CorrelatorInfo(avop1.operator_info,avop2.operator_info)
                    corr_name = repr(corr).replace(" ","-")
                    plh.add_correlator_subsection(repr(corr),self.proj_file_handler.corr_plot_file( corr_name, "pdf"),
                                                self.proj_file_handler.effen_plot_file( corr_name, "pdf"), index)

            self.proj_file_handler.remove_summary_files()
            if self.other_params['separate_mom']:
                loglevel = logging.getLogger().getEffectiveLevel()
                logging.getLogger().setLevel(logging.WARNING)
                pqdm([[self.proj_file_handler.summary_file(psq), i] for i,psq in enumerate(self.moms)], plh.compile_pdf, n_jobs=self.project_handler.nodes, argument_type="args")
                logging.getLogger().setLevel(loglevel)
                logging.info(f"Summary file saved to {self.proj_file_handler.summary_file('*')}.pdf.")
            else:
                plh.compile_pdf(self.proj_file_handler.summary_file()) 
                logging.info(f"Summary file saved to {self.proj_file_handler.summary_file()}.pdf.")
         
    #use class data to generate the channel plots
    def write_channel_plots_data(self, channel, plh):
        sigmond_util.write_channel_plots(self.averaged_operators[channel], plh, self.other_params['create_pickles'],
                            self.other_params['create_pdfs'] or self.other_params['create_summary'],self.proj_file_handler,
                            self.data[channel])
        
    #use estimates file to generate channel plots
    def write_channel_plots(self, channel, plh):
        sigmond_util.write_channel_plots(self.averaged_operators[channel], plh, self.other_params['create_pickles'],
                            self.other_params['create_pdfs'] or self.other_params['create_summary'],self.proj_file_handler)


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

    return operator_lib.Operator(averaged_channel.getGIOperator(obs_name, obs_id))


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