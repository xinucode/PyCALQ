import logging
import os
import pandas as pd
import numpy as np
import yaml

import general.task_manager as tm
import general.plotting_handler as ph
import general.particles
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
      rebin: 1                      #not required
      run_tag: ''                   #not required #default: ''
      sampling_mode: J              #not required
  - compare_files: []               #not required #default []
  - compare_rebin:                  #not required
      rebin_values: []              #required 
      run_tag: ''                   #not required #default: ''
      sampling_mode: J              #not required
      pivot_type: 0               #not required #default: 0
      t0: 8                         #required 
      tN: 5                         #required 
      tD: 18                        #required 
  - compare_tags:
      filetags: []                  #required
      sampling_mode: J              #not required
      pivot_type: 0               #not required #default: 0
      t0: 8                         #required 
      tN: 5                         #required 
      tD: 18                        #required 
      rebin: 1                      #not required
  figheight: 8                      #not required #default: 8
  figwidth: 15                      #not required #default: 15
  plot: true                        #required
  plot_deltaE: true                 #not required #default: True
  reference_particle: P             #not required #default: None
'''

#basic information about the rest mass particles (momentum is assumed zero)

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
            'max_level': 1000,
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
        
        #collect plot config infos
        self.compare_plots = []
        for compare_plot in task_configs['compare_plots']:
            root = list(compare_plot.keys())[0]
            plot_configs = compare_plot[root]
            sigmond_util.update_params(default_plot_configs,plot_configs)

            plot = {}
            #user defines the files #doesnt work
            if root=='compare_files':
                self.compare_plots.append(plot_configs)

                #compare spectrums of varying rebin values. If reference_particle is defined,
                # will plot relative error and chisqr of reference particle.
            elif root=='compare_rebin':
                for rebin in plot_configs['rebin_values']:
                    dataset_key = rf"$N_{{\textup{{bin}}}}={rebin}$"
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
                    sampling_mode = plot_configs['sampling_mode'] #+"-samplings"
                    key = proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].filekey(None, rebin, sampling_mode+"-samplings", rotate_type, tN, t0, tD, file_tag)
                    file = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].estimates_file(key)
                    if os.path.isfile(file):
                        plot[dataset_key] = file
                    else:
                        file2 = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].samplings_file( False, None, None, rebin, sampling_mode, 
                                                                                                    rotate_type, tN, t0, tD, 'levels'+file_tag)
                        if os.path.isfile(file2):
                            plot[dataset_key] = file2
                        else:
                            logging.warning(f"Could not find either '{file}' or '{file2}' for Nbin={rebin}.")

            #compare different pivots
            elif root=='compare_gevp':
                for pivot_set in plot_configs['gevp_values']:
                    tN = pivot_set['tN']
                    t0 = pivot_set['t0']
                    tD = pivot_set['tD']
                    rotate_type = 'SP'
                    if 'pivot_type' in pivot_set:
                        if pivot_set['pivot_type']:
                            rotate_type = 'RP'
                    dataset_key = f"{rotate_type}({tN},{t0},{tD})"
                    file_tag=''
                    if plot_configs['run_tag']:
                        file_tag='-'+plot_configs['run_tag']
                    sampling_mode = plot_configs['sampling_mode'] #+"-samplings"
                    rebin = plot_configs['rebin']
                    key = proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].filekey(None, rebin, sampling_mode+"-samplings", rotate_type, tN, t0, tD, file_tag)
                    file = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].estimates_file(key)
                    if os.path.isfile(file):
                        plot[dataset_key] = file
                    else:
                        file2 = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].samplings_file( False, None, None, rebin, sampling_mode, 
                                                                                                    rotate_type, tN, t0, tD, 'levels'+file_tag)
                        if os.path.isfile(file2):
                            plot[dataset_key] = file2
                        else:
                            logging.warning(f"Could not find either '{file}' or '{file2}' for pivot_set={dataset_key}.")

            #compare spectrums with different user defined tags. For those unexpected comparison, user-defined 
                #filetags are allowed in the spectrum task and then different tags can be compared here
            elif root=='compare_tags':
                for file_tag in plot_configs['filetags']:
                    dataset_key = file_tag
                    tN = plot_configs['tN']
                    t0 = plot_configs['t0']
                    tD = plot_configs['tD']
                    rotate_type = 'SP'
                    if 'pivot_type' in plot_configs:
                        if plot_configs['pivot_type']:
                            rotate_type = 'RP'
                    sampling_mode = plot_configs['sampling_mode']
                    # sampling_mode = sampling_mode+"-samplings"
                    rebin = plot_configs['rebin']
                    key = proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].filekey(None, rebin, sampling_mode+"-samplings", rotate_type, tN, t0, tD, '-'+file_tag)
                    file = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].estimates_file(key)
                    if os.path.isfile(file):
                        plot[dataset_key] = file
                    else:
                        file2 = self.proj_files_handler.all_tasks[tm.Task.fit_spectrum.name].samplings_file( False, None, None, rebin, sampling_mode, 
                                                                                                    rotate_type, tN, t0, tD, 'levels-'+file_tag)
                        if os.path.isfile(file2):
                            plot[dataset_key] = file2
                        else:
                            logging.warning(f"Could not find either '{file}' or '{file2}' for run_tag={file_tag}.")

            if plot:
                self.compare_plots.append(plot)
            else:
                logging.warning(f"Could not generate {root} plot.")
        
        #make yaml output
        logging.info(f"Full input written to '{os.path.join(proj_files_handler.log_dir(), 'full_input.yml')}'.")
        with open( os.path.join(proj_files_handler.log_dir(), 'full_input.yml'), 'w+') as log_file:
            yaml.dump({"general":general_configs, task_name: task_configs}, log_file)

    #this task is only designed to plot and compare
    def run( self ):
        pass

    
    def plot( self ):
        plh = ph.PlottingHandler()
        plh.create_fig(self.other_params['figwidth'], self.other_params['figheight'])
        plh.create_summary_doc("Compare Spectrums")

        datasets = {}
        #make the plots
        for iplot,plot in enumerate(self.compare_plots):
            #collect particles/channels involved
            particles = []
            for dataset in plot:
                if plot[dataset].endswith(".csv"):
                    if plot[dataset] not in datasets:
                        datasets[plot[dataset]] = pd.read_csv(plot[dataset])

                    #get channels
                    for i, row in datasets[plot[dataset]].iterrows():
                        particle = (row['isospin'], row['strangeness'])
                        particles.append(particle)

            particles = list(set(particles))

            #plot each channels on own graph, begin with summary plots
            plh.append_section(f"Summaries {iplot}")
            for particle in particles:
                plh.append_subsection(f"{particle[0]} S={particle[1]}")

                # set up error analysis data (for rebin analysis)
                study_particle=False
                error_analysis={}
                if self.other_params['reference_particle']:
                    if particle[0]==general.particles.data[self.other_params['reference_particle']]["isospin"]:
                        if particle[1]==general.particles.data[self.other_params['reference_particle']]["strangeness"]:
                            study_particle=True
                            error_analysis["x"] = []
                            error_analysis["val"] = []
                            error_analysis["err"] = []
                            error_analysis["chisqrdof"] = []

                #add each spectrum set one at a time
                plh.clf()
                nothing = True
                this_energy_key = self.energy_key
                for id,dataset in enumerate(plot):
                    df = datasets[plot[dataset]]
                    if f'{self.energy_key} value' not in df:
                        this_energy_key = 'ecm'
                for id,dataset in enumerate(plot):

                    #set up spectrum data
                    indexes = []
                    levels = []
                    errs = []
                    irreps = []
                    df = datasets[plot[dataset]]

                    for i, row in df[(df['isospin']==particle[0])&(df['strangeness']==particle[1])].iterrows():
                        irrep = (row['irrep'],row['momentum'])
                        if irrep not in irreps:
                            irreps.append(irrep)

                    levels_key = 'fit level'
                    if type(df[levels_key][0])==np.float64:
                        levels_key = 'rotate level'


                    #sort irreps
                    irreps.sort(key=lambda x: x[1]+psettings.alphabetical[x[0]])

                    #collect spectrum and error analysis data
                    for i, row in df[(df['isospin']==particle[0])&(df['strangeness']==particle[1])].iterrows():
                        if not np.isnan(row[f'{this_energy_key} value']) and row[levels_key]<=self.other_params['max_level']:
                            irrep = (row['irrep'],row['momentum'])
                            indexes.append(irreps.index(irrep))
                            levels.append(row[f'{this_energy_key} value'])
                            errs.append(row[f'{this_energy_key} error'])
                            nothing = False
                            if study_particle and row['momentum']==0 and row[levels_key]==0:
                                error_analysis['x'].append(dataset)
                                error_analysis['val'].append(row[f'{this_energy_key} value'])
                                error_analysis['err'].append(row[f'{this_energy_key} error'])
                                error_analysis['chisqrdof'].append(row['chisqrdof'])

                    #plot spectrum dataset
                    if levels:
                        if this_energy_key=="ecm":
                            plh.summary_plot(indexes,levels,errs,irreps, None, [], dataset, id, len(plot))
                        else:
                            plh.summary_plot(indexes,levels,errs,irreps, self.other_params['reference_particle'], [], dataset, id, len(plot))

                #finalize spectrum comparisons
                if not nothing:
                    strangeness = particle[1]
                    if particle[1]<0:
                        strangeness = -particle[1]
                        strangeness = f"m{particle[1]}"
                    filekey = f"{self.compare_plots.index(plot)}-{particle[0]}_S{strangeness}"
                    plh.save_pdf(self.proj_files_handler.summary_plot_file("pdf",filekey))
                    logging.info(f"Comparison plot saved to '{self.proj_files_handler.summary_plot_file('pdf',filekey)}'.")
                    plh.add_single_plot(self.proj_files_handler.summary_plot_file("pdf",filekey))
                
                #generate error analysis plot
                if study_particle and error_analysis['x']:
                    plh.clf()
                    relerrs = np.array(error_analysis['err'])/np.array(error_analysis['val'])
                    relerrs /= relerrs[0]
                    plh.show_trend(error_analysis['x'],relerrs,"$R_N/R_1$")
                    plh.show_trend(error_analysis['x'],error_analysis['chisqrdof'],r"$\chi^2/\textup{dof}$", True)
                    plh.save_pdf(os.path.join(self.proj_files_handler.plot_dir("pdfs"),f"error_analysis-{self.compare_plots.index(plot)}.pdf"))
                    plh.add_single_plot(os.path.join(self.proj_files_handler.plot_dir("pdfs"),f"error_analysis-{self.compare_plots.index(plot)}.pdf"))

            #plot shifts if asked for and given
            if self.other_params['plot_deltaE']:
                plh.append_section(f"Energy shifts {iplot}")
                for particle in particles:
                    plh.append_subsection(f"{particle[0]} S={particle[1]}")
                    irreps = {}
                    for id,dataset in enumerate(plot):
                        df = datasets[plot[dataset]]
                        levels_key = 'fit level'
                        if type(df[levels_key][0])==np.float64:
                            levels_key = 'rotate level'
                        if 'dElab value' not in df:
                            continue
                        for i, row in df[(df['isospin']==particle[0])&(df['strangeness']==particle[1])].iterrows():
                            if row['momentum'] not in irreps:
                                irreps[row['momentum']] = []

                            if row[levels_key]>self.other_params['max_level']:
                                continue

                            # print(dataset,row['momentum'],irreps[row['momentum']],row['irrep'],row[levels_key])
                            if row['irrep'] not in irreps[row['momentum']]:
                                irreps[row['momentum']].append(row['irrep'])
                            elif irreps[row['momentum']].index(row['irrep'])+row[levels_key]>=len(irreps[row['momentum']]):
                                irreps[row['momentum']].append(row['irrep'])
                            elif irreps[row['momentum']][irreps[row['momentum']].index(row['irrep'])+int(row[levels_key])]!=row['irrep']:
                                irreps[row['momentum']].insert(irreps[row['momentum']].index(row['irrep']),row['irrep'])


                        #sort irreps
                        [irreps[irrep].sort(key=lambda x: psettings.alphabetical[x]) for irrep in irreps]

                    #each momentum is on a different plot
                    for mom in irreps:
                        plh.clf()
                        nothing = True
                        for id,dataset in enumerate(plot):
                            df = datasets[plot[dataset]]
                            levels_key = 'fit level'
                            if type(df[levels_key][0])==np.float64:
                                levels_key = 'rotate level'
                            indexes = []
                            levels = []
                            errs = []
                            split_irreps = []
                            if 'dElab value' not in df:
                                continue
                            for i, row in df[(df['isospin']==particle[0])&(df['strangeness']==particle[1])&(df['momentum']==mom)].iterrows():
                                if not np.isnan(row[f'dElab value']) and row[levels_key]<=self.other_params['max_level']:
                                    irrep = row['irrep']
                                    fit_level = row[levels_key]
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
                    