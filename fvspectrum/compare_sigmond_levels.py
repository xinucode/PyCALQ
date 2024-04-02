import logging
import os
import pandas as pd
import numpy as np
import yaml

import general.task_manager as tm
import general.plotting_handler as ph
import fvspectrum.sigmond_util as sigmond_util
import fvspectrum.spectrum_plotting_settings.settings as psettings

doc = '''
compare_spectrums - using the fit results of the fit_spectrum task, plots
    the spectrums side by side for comparisons and puts in summary document. WIP.

task input
---
compare_spectrums:              #required
  compare_plots:                #required # list of comparison plot types
  - compare_gevp:               #not required
      gevp_values:              #required for "compare_gevp", list of pivot configs
      - t0: 8                       #required
        tD: 16                      #required
        tN: 5                       #required
        pivot_type: 0               #not required #default: 0
      - t0: 8
        tD: 18
        tN: 5
      ...
      rebin: 1                      #required
      run_tag: ''                   #not required #default: ''
      sampling_mode: J              #required
  figheight: 8                      #not required #default: 8
  figwidth: 15                      #not required #default: 15
  plot: true                        #required
  plot_deltaE: true                 #not required #default: True
  reference_particle: P             #not required #default: None
'''

class CompareLevels:
    @property
    def info(self):
        return doc

    def __init__( self, task_name, proj_files_handler, general_configs, task_configs ):
        self.proj_files_handler= proj_files_handler
        self.other_params = {
            # 'create_pdfs': True,
            # 'create_pickles': True,
            # 'create_summary': True,
            'plot': True,
            'figwidth':15,
            'figheight':8,
            'reference_particle': None,
            'plot_deltaE': True,
        }
        sigmond_util.update_params(self.other_params,task_configs)
        if not self.other_params['plot']:
            logging.warning(f"No plots requested, task {task_name} does nothing.")
            return
        
        self.energy_key = 'ecm'
        if self.other_params['reference_particle']!=None:
            self.energy_key = 'ecm_ref'

        
        default_plot_configs = {
            "rebin": 1,
            "sampling_mode": 'J',
            "run_tag": "",
        }
        
        self.compare_plots = []
        for compare_plot in task_configs['compare_plots']:
            root = list(compare_plot.keys())[0]
            plot_configs = compare_plot[root]
            sigmond_util.update_params(default_plot_configs,plot_configs)
            if root=='compare_files':
                self.compare_plots.append(plot_configs)
            if root=='compare_rebin':
                plot = {}
                for rebin in plot_configs['rebin_values']:
                    dataset_key = rf"$N_{{\textup{{bin}}}}={rebin}"
                    file_tag=''
                    if plot_configs['run_tag']:
                        file_tag='-'+plot_configs['run_tag']
                    tN = plot_configs['tN']
                    t0 = plot_configs['t0']
                    tD = plot_configs['tD']
                    rotate_type = 'SP'
                    if 'pivot_type' in plot_configs:
                        if plot_configs['pivot_type']:
                            rotate_type = 'RP'
                    sampling_mode = plot_configs['sampling_mode']
                    key = proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].filekey(None, rebin, sampling_mode, rotate_type, tN, t0, tD, file_tag)
                    file = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].estimates_file(key)
                    if os.path.isfile(file):
                        plot[dataset_key] = file
                    else:
                        file2 = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].samplings_file( False, None, None, rebin, sampling_mode, 
                                                                                                    rotate_type, tN, t0, tD, file_tag)
                        if os.path.isfile(file):
                            plot[dataset_key] = file2
                        else:
                            logging.warning(f"Could not find either '{file}' or '{file2}' for Nbin={rebin}.")
                if plot:
                    self.compare_plots.append(plot)
                else:
                    logging.warning(f"Could not generate rebin comparison for rebin values: {plot_configs['rebin_values']}.")
            if root=='compare_gevp':
                plot = {}
                for pivot_set in plot_configs['gevp_values']:
                    tN = pivot_set['tN']
                    t0 = pivot_set['t0']
                    tD = pivot_set['tD']
                    rotate_type = 'SP'
                    if 'pivot_type' in pivot_set:
                        if pivot_set['pivot_type']:
                            rotate_type = 'RP'
                    dataset_key = f"({tN},{t0},{tD})"
                    file_tag='' #???
                    if plot_configs['run_tag']:
                        file_tag='-'+plot_configs['run_tag']
                    sampling_mode = plot_configs['sampling_mode']+"-samplings"
                    rebin = plot_configs['rebin']
                    key = proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].filekey(None, rebin, sampling_mode, rotate_type, tN, t0, tD, file_tag)
                    file = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].estimates_file(key)
                    if os.path.isfile(file):
                        plot[dataset_key] = file
                    else:
                        file2 = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].samplings_file( False, None, None, rebin, sampling_mode, 
                                                                                                    rotate_type, tN, t0, tD, file_tag)
                        if os.path.isfile(file):
                            plot[dataset_key] = file2
                        else:
                            logging.warning(f"Could not find either '{file}' or '{file2}' for Nbin={rebin}.")
                if plot:
                    self.compare_plots.append(plot)
                else:
                    logging.warning(f"Could not generate gevp comparison for rebin values: {plot_configs['gevp_values']}.")
        
        #make yaml output
        logging.info(f"Full input written to '{os.path.join(proj_files_handler.log_dir(), 'full_input.yml')}'.")
        with open( os.path.join(proj_files_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
            yaml.dump({"general":general_configs, task_name: task_configs}, log_file)

    def run( self ):
        pass

    def plot( self ):
        plh = ph.PlottingHandler()
        plh.create_fig(self.other_params['figwidth'], self.other_params['figheight'])
        plh.create_summary_doc("Compare Spectrums")

        # self.compare_plots = list(set(self.compare_plots))
        datasets = {}
        for plot in self.compare_plots:
            for dataset in plot:
                if plot[dataset].endswith(".csv"):
                    if plot[dataset] not in datasets:
                        datasets[plot[dataset]] = pd.read_csv(plot[dataset])
                    particles = []

                    #get channels
                    for i, row in datasets[plot[dataset]].iterrows():
                        particle = (row['isospin'], row['strangeness'])
                        particles.append(particle)

            
            particles = list(set(particles))

            #plot each channels on own graph
            plh.append_section("Summaries")
            for particle in particles:
                plh.append_subsection(f"{particle[0]} S={particle[1]}")
                plh.clf()
                nothing = True
                for id,dataset in enumerate(plot):
                    indexes = []
                    levels = []
                    errs = []
                    irreps = []
                    df = datasets[plot[dataset]]
                    for i, row in df[(df['isospin']==particle[0])&(df['strangeness']==particle[1])].iterrows():
                        irrep = (row['irrep'],row['momentum'])
                        if irrep not in irreps:
                            irreps.append(irrep)

                    #sort irreps
                    irreps.sort(key=lambda x: x[1]+psettings.alphabetical[x[0]])

                    for i, row in df[(df['isospin']==particle[0])&(df['strangeness']==particle[1])].iterrows():
                        if not np.isnan(row[f'{self.energy_key} value']):
                            irrep = (row['irrep'],row['momentum'])
                            indexes.append(irreps.index(irrep))
                            levels.append(row[f'{self.energy_key} value'])
                            errs.append(row[f'{self.energy_key} error'])
                            nothing = False
                    if levels:
                        plh.summary_plot(indexes,levels,errs,irreps, self.other_params['reference_particle'], [], dataset, id, len(plot))
                if not nothing:
                    strangeness = particle[1]
                    if particle[1]<0:
                        strangeness = -particle[1]
                        strangeness = f"m{particle[1]}"
                    filekey = f"{self.compare_plots.index(plot)}-{particle[0]}_S{strangeness}"
                    plh.save_pdf(self.proj_files_handler.summary_plot_file("pdf",filekey))
                    logging.info(f"Comparison plot saved to '{self.proj_files_handler.summary_plot_file('pdf',filekey)}'.")
                    plh.add_single_plot(self.proj_files_handler.summary_plot_file("pdf",filekey))

            if self.other_params['plot_deltaE']:
                plh.append_section("Energy shifts")
                for particle in particles:
                    plh.append_subsection(f"{particle[0]} S={particle[1]}")
                    irreps = {}
                    for id,dataset in enumerate(plot):
                        df = datasets[plot[dataset]]
                        for i, row in df[(df['isospin']==particle[0])&(df['strangeness']==particle[1])].iterrows():
                            if row['momentum'] not in irreps:
                                irreps[row['momentum']] = []
                            if row['irrep'] not in irreps[row['momentum']]:
                                irreps[row['momentum']].append(row['irrep'])
                            elif irreps[row['momentum']].index(row['irrep'])+row['fit level']>=len(irreps[row['momentum']]):
                                irreps[row['momentum']].append(row['irrep'])
                            elif irreps[row['momentum']][irreps[row['momentum']].index(row['irrep'])+row['fit level']]!=row['irrep']:
                                irreps[row['momentum']].insert(irreps[row['momentum']].index(row['irrep']),row['irrep'])


                        #sort irreps
                        [irreps[irrep].sort(key=lambda x: psettings.alphabetical[x]) for irrep in irreps]

                    for mom in irreps:
                        plh.clf()
                        nothing = True
                        for id,dataset in enumerate(plot):
                            df = datasets[plot[dataset]]
                            indexes = []
                            levels = []
                            errs = []
                            split_irreps = []
                            for i, row in df[(df['isospin']==particle[0])&(df['strangeness']==particle[1])&(df['momentum']==mom)].iterrows():
                                if not np.isnan(row[f'dElab value']):
                                    irrep = row['irrep']
                                    fit_level = row['fit level']
                                    split_irreps.append((irrep, mom, fit_level))
                                    indexes.append(irreps[mom].index(irrep)+fit_level)
                                    levels.append(row[f'dElab value'])
                                    errs.append(row[f'dElab error'])
                                    nothing = False
                            if levels:
                                plh.summary_plot(indexes,levels,errs,split_irreps, None, [], dataset, id, len(plot), True)
                        if not nothing:
                            strangeness = particle[1]
                            if particle[1]<0:
                                strangeness = -particle[1]
                                strangeness = f"m{particle[1]}"
                            filekey = f"{self.compare_plots.index(plot)}-{particle[0]}_S{strangeness}_PSQ{mom}"
                            plh.save_pdf(self.proj_files_handler.summary_plot_file("pdf",filekey))
                            logging.info(f"Comparison plot saved to '{self.proj_files_handler.summary_plot_file('pdf',filekey)}'.")
                            plh.add_single_plot(self.proj_files_handler.summary_plot_file("pdf",filekey))

        plh.compile_pdf(self.proj_files_handler.summary_file())
        logging.info(f"Summary file saved to {self.proj_files_handler.summary_file()}.pdf")
                    