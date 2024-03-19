import logging
import os
import yaml
import pandas as pd

import sigmond
import fvspectrum.sigmond_util as sigmond_util
import general.plotting_handler as ph
# import fvspectrum.sigmond_info.sigmond_input as sigmond_input
# import fvspectrum.sigmond_info.sigmond_info as sigmond_info
# import fvspectrum.sigmond_operator_info.operator as operator
# from fvspectrum.sigmond_data_handling.correlator_data import CorrelatorData
# import fvspectrum.sigmond_data_handling.data_handler as data_handler
from sigmond_scripts.analysis.sigmond_info import sigmond_info, sigmond_input
from sigmond_scripts.analysis.operator_info import operator
from sigmond_scripts.analysis.data_handling import data_handler
from sigmond_scripts.analysis.data_handling.correlator_data import CorrelatorData

doc = '''
general:
  ensemble_id: cls21_c103               #required
  ensembles_file: /home/sarahski/latticeQCD/pycalq/pycalq.new_sigmond_task/fvspectrum/sigmond_utils/ensembles.xml #generated by PyCALQ
  project_dir: /latticeQCD/raid3/sarahski/lqcd/C103_R005/test_pycalq_project #required
  sampling_info:                        #not required 
    mode: Jackknife                     #default Jackknife
  tweak_ensemble:                       #not required
    omissions: []                       #default []
    rebin: 1                            #default 1
rotate_corrs:                           #required
  t0: 5                                 #required
  tD: 10                                #required
  tN: 5                                 #required
  averaged_input_correlators_dir: {project_dir}/1average_corrs/data/bins #not required #default {project_dir}/1average_corrs/data/bins 
  create_pdfs: true                     #not required #default true
  create_pickles: true                  #not required #default true
  create_summary: true                  #not required #default true
  erase_original_matrix_from_memory: true #not required #default true
  figheight: 6                          #not required #default 6
  figwidth: 8                           #not required #default 8
  generate_estimates: true              #not required #default true
  ignore_missing_correlators: true      #not required #default true
  max_condition_number: 50              #not required #default 50
  only:                                 #not required
  - psq=0
  - isosinglet S=0 E PSQ=3
  omit:                                 #not required (overridden by 'only' setting)
  - psq=0
  - isosinglet S=0 E PSQ=3
  pivot_type: 0                         #not required #default 0; 0 - single pivot, 1 - rolling pivot
  plot: true                            #not required #default true
  precompute: true                      #not required #default true
  rotate_by_samplings: true             #not required #default true; otherwise rotate by bins
  tmax: 25                              #not required #default 25
  tmin: 2                               #not required #default 2
  used_averaged_bins: true              #not required #default true
'''
#to do: add specilized pivots for individual irreps, get average directory without hardcode
    #address tmin/tmax problem
class SigmondRotateCorrs:
    @property
    def info(self):
        return doc

    def rotated_corrs_dir( self, binned ): #add rotation info, and then average info
        if binned:
            subdir = 'bins'
        else:
            subdir = 'samples'
        return self.proj_dir_handler.data_dir(subdir)
        
    def rotated_corrs_file( self, binned, channel = None ): #add rotation info, and then average info
        if channel:
            return os.path.join(self.rotated_corrs_dir(binned),f"rotated_corrs.hdf5[{channel}]")
        else:
            return os.path.join(self.rotated_corrs_dir(binned),f"rotated_corrs.hdf5")

    def pivot_file( self, channel=None ):
        if channel:
            return os.path.join(self.proj_dir_handler.data_dir('pivots'),f"pivot_information.hdf5[{channel}]")
        else:
            return os.path.join(self.proj_dir_handler.data_dir('pivots'),f"pivot_information.hdf5")
        
    def rotated_correstimates_file(self, corr):
        return os.path.join(self.proj_dir_handler.data_dir('estimates'), f"{corr}_correlator_estimates.csv")
    def rotated_effestimates_file(self, corr):
        return os.path.join(self.proj_dir_handler.data_dir('estimates'), f"{corr}_effenergy_estimates.csv")
    
    def corr_plot_file(self,corr, ptype):
        return os.path.join(self.proj_dir_handler.plot_dir(f"{ptype}s"), f"{corr}_correlator.{ptype}")

    def effen_plot_file(self,corr, ptype):
        return os.path.join(self.proj_dir_handler.plot_dir(f"{ptype}s"), f"{corr}_effenergy.{ptype}")
    
    @property
    def summary_file(self):
        return os.path.join(self.proj_dir_handler.plot_dir(), f"{self.task_name}_summary") #add channel? project name?
        
    def __init__( self, task_name, proj_dir_handler, general_configs, task_configs ):
        self.task_name = task_name
        self.proj_dir_handler= proj_dir_handler

        self.project_info = sigmond_util.setup_project(general_configs)

        #other params
        self.other_params = {
            'create_pdfs': True,
            'create_pickles': True,
            'create_summary': True,
            'plot': True,
            'figwidth':8,
            'figheight':6,
            'tmin':2,
            'tmax':25,
            'erase_original_matrix_from_memory': True,
            'ignore_missing_correlators': True,
            'generate_estimates': True,
            'rotate_by_samplings': True,
            'pivot_type': 0, #0 - single; 1 - rolling
            'precompute': True,
            'max_condition_number': 50,
            'used_averaged_bins': True #otherwise samplings
        }
        sigmond_util.update_params(self.other_params,task_configs) #update other_params with task_params, 
                                                                        #otherwise fill in missing task params

        #get averaged channels 
        this_data_handler = data_handler.DataHandler(self.project_info)
        averaged_directories = []
        averaged_bin_dir = os.path.join(self.project_info.project_dir,"1average_corrs","data",'bins')
        averaged_samples_dir = os.path.join(self.project_info.project_dir,"1average_corrs","data",'samples')
        if 'averaged_input_correlators_dir' in task_configs:
            if type(task_configs['averaged_input_correlators_dir'])==list:
                averaged_directories+=task_configs['averaged_input_correlators_dir']
            else:
                averaged_directories.append(task_configs['averaged_input_correlators_dir'])
        elif self.other_params['used_averaged_bins'] and os.path.isdir(averaged_bin_dir):
            averaged_directories.append(os.path.join("1average_corrs","data",'bins'))
            task_configs['averaged_input_correlators_dir'] = averaged_bin_dir
        elif os.path.isdir(averaged_samples_dir):
            averaged_directories.append(os.path.join("1average_corrs","data",'samples'))
            task_configs['averaged_input_correlators_dir'] = averaged_samples_dir

        this_data_handler.rel_averaged_datadir = averaged_directories
        self.averaged_data_files = this_data_handler.findAveragedData()

        #establish rotation parameters
        for param in ['t0','tN','tD']:
            if param not in task_configs:
                logging.critical(f"No default '{param}' set for '{task_name}'. Ending task.")
        
        self.t0 = task_configs['t0']
        self.tN = task_configs['tN']
        self.tD = task_configs['tD']

        #do only channels specified
        self.data_handler = this_data_handler
        self.channels = self.data_handler.averaged_channels[:]
        if 'only' in task_configs:
            only_moms = []
            only_channels = []
            for item in task_configs['only']:
                if item.startswith('PSQ='):
                    only_moms.append(int(item.replace('PSQ=',"")))
                elif item.startswith('psq='):
                    only_moms.append(int(item.replace('psq=',"")))
                else:
                    only_channels.append(item)
            final_channels = []
            for channel in self.channels:
                if channel.psq in only_moms or str(channel) in only_channels:
                    final_channels.append(channel)
                else:
                    logging.info(f'Channel {str(channel)} omitted due to "only" setting.')
            self.channels = final_channels
        elif 'omit' in task_configs:
            omit_moms = []
            omit_channels = []
            for item in task_configs['omit']:
                if item.startswith('PSQ='):
                    omit_moms.append(int(item.replace('PSQ=',"")))
                elif item.startswith('psq='):
                    omit_moms.append(int(item.replace('psq=',"")))
                else:
                    omit_channels.append(item)
            final_channels = []
            for channel in self.channels:
                if channel.psq in omit_moms or str(channel) in omit_channels:
                    logging.info(f'Channel {str(channel)} omitted due to "omit" setting.')
                else:
                    final_channels.append(channel)
            self.channels = final_channels

        #remove unqualified channels
        rm_channels = []
        for channel in self.channels:
            if len(self.data_handler.getChannelOperators(channel))<2:
                rm_channels.append(channel)
                logging.info(f"Skipping {str(channel)} because there is an insufficient number of operators.")
        for channel in rm_channels:
            self.channels.remove(channel)

        task_configs['rotated_channels'] = []
        for channel in self.channels:
            task_configs['rotated_channels'].append(str(channel))

        #these wont change, all correlators can and should be considered hermetian, 
        #       and time separation is a cosmetic parameter
        self.hermitian = True
        self.time_separation = 1

        #these can and will matter but only for special cases. Will need extra care when coding up. 
        #only coding up if we come across an instance of needing such
        self.subtract_vev = False
        self.vev_const = 0.0
        self.effective_energy_type = 0 #0=TimeForward, 1=TimeSymmetric, 2=TimeBackward?
        
        if self.other_params['pivot_type']>1 or self.other_params['pivot_type']<0:
            logging.critical(f"Parameter 'pivot_type' set to unknown pivot type, select 0 for single pivot, and 1 for rolling pivot.")

        if not self.other_params['create_pdfs'] and not self.other_params['create_pickles'] and not self.other_params['create_summary']:
            self.other_params['plot'] = False
        
        #make yaml output
        logging.info(f"Full input written to '{os.path.join(proj_dir_handler.log_dir(), 'full_input.yml')}'.")
        with open( os.path.join(proj_dir_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
            yaml.dump({"general":general_configs, task_name: task_configs}, log_file)

    def run( self ):
        #split for parallel jobs
        task_input = sigmond_input.SigmondInput(
            os.path.basename(self.project_info.project_dir),
            self.project_info.bins_info,
            self.project_info.sampling_info,
            self.project_info.ensembles_file,
            self.averaged_data_files,
            "temp1.xml",
            os.path.join(self.proj_dir_handler.log_dir(),"sigmond_rotation_log.xml"), #actually creating this
            self.other_params['precompute'],
            None,
        )
        if self.other_params['pivot_type']==0:
            pivot_string = 'single_pivot'
        elif self.other_params['pivot_type']==1:
            pivot_string = 'rolling_pivot'

        file_created = False
        self.nlevels = {}
        skip_channels = []
        for channel in self.channels:
            if file_created:
                wmode = sigmond.WriteMode.Update
                overwrite = False
            else:
                wmode = sigmond.WriteMode.Overwrite
                overwrite = True
            operators = set([op.operator_info for op in self.data_handler.getChannelOperators(channel)])
            self.nlevels[channel] = len(operators)
            if len(operators)<2:
                logging.info(f"Skipping {str(channel)} because there is an insufficient number of operators.")
                skip_channels.append(channel)
                continue

            if self.other_params['rotate_by_samplings']:
                rotate_mode = "samplings_all"
            else:
                rotate_mode = "bins"
            task_input.doCorrMatrixRotation(
                sigmond_info.PivotInfo(pivot_string,norm_time=self.tN,metric_time=self.t0,
                                    diagonalize_time=self.tD,max_condition_number=self.other_params['max_condition_number']),
                sigmond_info.RotateMode(rotate_mode),
                sigmond.CorrelatorMatrixInfo(operators,self.hermitian,self.subtract_vev),
                operator.Operator( channel.getRotatedOp() ),
                self.other_params['tmin'],
                self.other_params['tmax'],
                rotated_corrs_filename=self.rotated_corrs_file(not self.other_params['rotate_by_samplings'],repr(channel)),
                file_mode=wmode,
                pivot_filename=self.pivot_file(repr(channel)),
                pivot_overwrite=overwrite,
            )
            file_created = True

        for channel in skip_channels:
            self.channels.remove(channel)

        setuptaskhandler = sigmond.XMLHandler()
        setuptaskhandler.set_from_string(task_input.to_str())

        with open(os.path.join(self.proj_dir_handler.log_dir(),'sigmond_rotation_input.xml'), 'w+') as f:
            f.write(setuptaskhandler.output(1))
        logging.info(f"Sigmond input file written to {os.path.join(self.proj_dir_handler.log_dir(),'sigmond_rotation_input.xml')}")

        taskhandler = sigmond.TaskHandler(setuptaskhandler)
        taskhandler.do_batch_tasks(setuptaskhandler)
        del taskhandler

        if os.path.isfile(self.pivot_file()):
            logging.info(f"Pivot matrix written to {self.pivot_file()}.") #add FAIR
        if os.path.isfile(self.rotated_corrs_file(not self.other_params['rotate_by_samplings'])):
            logging.info(f"Rotated correlators written to {self.rotated_corrs_file(not self.other_params['rotate_by_samplings'])}.") #add FAIR
        logging.info(f"Log file written to {os.path.join(self.proj_dir_handler.log_dir(),'sigmond_rotation_log.xml')}")

        # data_files = CorrelatorData()
        # self.data_handler._averaged_data = data_files
        if self.other_params['plot'] or self.other_params['generate_estimates']:
            self.data_handler.rel_rotated_datadir = [self.rotated_corrs_dir(not self.other_params['rotate_by_samplings'])]
            self.data_handler.findRotatedData()
            mcobs_handler, mcobs_get_handler = sigmond_util.get_mcobs_handlers(self.data_handler, self.project_info)

        if self.other_params['plot'] and not self.other_params['generate_estimates']:
            self.rotated_estimates = {}
        if self.other_params['generate_estimates'] or self.other_params['plot']:
            logging.info(f"Generating estimates for {self.proj_dir_handler.data_dir('estimates')}...")
            for channel in self.channels:
                if self.other_params['plot'] and not self.other_params['generate_estimates']:
                    self.rotated_estimates[str(channel)] = {}
                logging.info(f"\tGenerating estimates for channel {str(channel)}...")
                for i in range(self.nlevels[channel]):
                    for j in range(self.nlevels[channel]):
                        rop1 = operator.Operator( channel.getRotatedOp(i) )
                        rop2 = operator.Operator( channel.getRotatedOp(j) )
                        corr = sigmond.CorrelatorInfo(rop1.operator_info,rop2.operator_info)
                        corr_name = repr(corr).replace(" ","-")
                        if self.other_params['plot'] and not self.other_params['generate_estimates']:
                            self.rotated_estimates[str(channel)][corr] = {}
                        
                        estimates = sigmond.getCorrelatorEstimates(mcobs_handler,corr,self.hermitian,self.subtract_vev,sigmond.ComplexArg.RealPart, 
                                                               self.project_info.sampling_info.getSamplingMode())
                        if self.other_params['generate_estimates']:
                            if estimates:
                                sigmond_util.estimates_to_csv(estimates, self.rotated_correstimates_file(corr_name) )
                            else:
                                if os.path.exists(self.rotated_correstimates_file(corr_name)):
                                    os.remove(self.rotated_correstimates_file(corr_name))
                        elif self.other_params['plot']:
                            self.rotated_estimates[str(channel)][corr]["corr"] = sigmond_util.estimates_to_df(estimates)
                        if not estimates:
                            logging.warning(f"No data found for {repr(corr)}.")

                        estimates = sigmond.getEffectiveEnergy(mcobs_handler,corr,self.hermitian,self.subtract_vev,sigmond.ComplexArg.RealPart, 
                                                           self.project_info.sampling_info.getSamplingMode(),self.time_separation,self.effective_energy_type,self.vev_const)
                        if self.other_params['generate_estimates']:
                            if estimates:
                                sigmond_util.estimates_to_csv(estimates, self.rotated_effestimates_file(corr_name) )
                            else:
                                if os.path.exists(self.rotated_effestimates_file(corr_name)):
                                    os.remove(self.rotated_effestimates_file(corr_name))
                        elif self.other_params['plot']:
                            self.rotated_estimates[str(channel)][corr]["effen"] = sigmond_util.estimates_to_df(estimates)

        
    def plot( self ):
        if self.other_params['plot']:
            logging.info(f"Saving plots to directory {self.proj_dir_handler.plot_dir()}...")
        else:
            logging.info(f"No plots requested.")
            return
        
        plh = ph.PlottingHandler()
        plh.create_fig(self.other_params['figwidth'], self.other_params['figheight'])
        if self.other_params['create_summary']:
            plh.create_summary_doc("Rotated Correlators")
        

        for channel in self.channels:
            logging.info(f"\tGenerating plots for channel {str(channel)}...")
            if self.other_params['create_summary']:
                plh.append_section(str(channel))
            corr_order = [(i,i) for i in range(self.nlevels[channel])]
            for i in range(self.nlevels[channel]):
                for j in range(self.nlevels[channel]):
                    if j is not i:
                        corr_order.append((i,j))
            for i,j in corr_order:
                rop1 = operator.Operator( channel.getRotatedOp(i) )
                rop2 = operator.Operator( channel.getRotatedOp(j) )
                corr = sigmond.CorrelatorInfo(rop1.operator_info,rop2.operator_info)
                corr_name = repr(corr).replace(" ","-")

                if not self.other_params['generate_estimates']:
                    df = self.rotated_estimates[str(channel)][corr]["corr"]
                    if df.empty:
                        continue
                else:
                    if os.path.exists(self.rotated_correstimates_file(corr_name)):
                        df = pd.read_csv(self.rotated_correstimates_file(corr_name))
                    else: 
                        continue
                
                plh.clf()
                plh.correlator_plot(df, 0, rop1, rop2) #0 for regular corr plot

                if self.other_params['create_pickles']:
                    plh.save_pickle(self.corr_plot_file( corr_name, "pickle"))
                if self.other_params['create_pdfs'] or self.other_params['create_summary']:
                    plh.save_pdf(self.corr_plot_file( corr_name, "pdf"))

                try:
                    if not self.other_params['generate_estimates']:
                        df = self.rotated_estimates[str(channel)][corr]["effen"]
                    else:
                        df = pd.read_csv(self.rotated_effestimates_file(corr_name))

                    plh.clf()
                    plh.correlator_plot(df, 1, rop1, rop2) #1 for effective energy plot

                    if self.other_params['create_pickles']:
                        plh.save_pickle(self.effen_plot_file( corr_name, "pickle"))
                    if self.other_params['create_pdfs'] or self.other_params['create_summary']:
                        plh.save_pdf( self.effen_plot_file(corr_name, "pdf")) 
                except FileNotFoundError as err:
                    pass
                except KeyError as err:
                    pass

                if self.other_params['create_summary']:
                    if i==j:
                        corr_title = str(rop1)
                    else:
                        corr_title = repr(corr)
                    plh.add_correlator_subsection(corr_title,self.corr_plot_file( corr_name, "pdf"),
                                                    self.effen_plot_file( corr_name, "pdf"))
                        
        if self.other_params['create_summary']:
            plh.compile_pdf(self.summary_file) 
            logging.info(f"Summary file saved to {self.summary_file}.")
        