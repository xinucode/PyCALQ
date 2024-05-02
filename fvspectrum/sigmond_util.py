import logging
import os
import xmltodict
import xml.etree.ElementTree as ET
import matplotlib
import shutil
from typing import NamedTuple
import pandas as pd
import numpy as np
import scipy as sp
import time
from threading import Thread

import sigmond
from sigmond_scripts import data_files, data_handler
from sigmond_scripts import fit_info, sigmond_info, sigmond_log
from sigmond_scripts import channel

#project info class drew designed for keeping all
    #important information for correlator analysis
class ProjectInfo(NamedTuple):
  project_dir: str
  raw_data_dirs: list
  ensembles_file: str
  echo_xml: bool
  bins_info: sigmond.MCBinsInfo
  sampling_info: sigmond.MCSamplingInfo
  data_files: data_files.DataFiles
  precompute: bool
  latex_compiler: str

#check that raw data files are real and not within project
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

#from the ensemble file, retrieve ensemble info indicated by 'ensemble_id'
    #general_params
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

#set up ProjectInfo class based on general parameters and raw data list
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

    datafiles = data_files.DataFiles()

    return ProjectInfo( project_dir=general_params['project_dir'], raw_data_dirs=raw_data_files, ensembles_file=ensemble_file_path,
            echo_xml=False, bins_info=bins_info, sampling_info=sampling_info, data_files=datafiles,
            precompute=True, latex_compiler=None)

#check if latex is available on the system
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

#for each default param in other_params, the user set parameters
    #in task_params will override the original value
def update_params( other_params, task_params):
    for param in other_params:
        if param in task_params:
            other_params[param] = task_params[param]
        else:
            task_params[param] = other_params[param]

#get all data handlers, data_handler for files and general channel info
    #mcobs_handler for samplings and bins themselves, must return mcobs_get_handler
    #else mcobs_handler does not work
def get_data_handlers(project_info):
    this_data_handler = data_handler.DataHandler(project_info)
    mcobs_handler, mcobs_get_handler = get_mcobs_handlers(this_data_handler, project_info)
    return this_data_handler, mcobs_handler, mcobs_get_handler

#get mcobs_handler for managing bin and sample data, must return mcobs_get_handler
    #else mcobs_handler does not work
def get_mcobs_handlers(this_data_handler, project_info, additional_sampling_files = []):
    mcobs_handler_init = ET.Element("MCObservables")
    if this_data_handler.raw_data_files.bl_corr_files or this_data_handler.averaged_data_files.bl_corr_files or this_data_handler.rotated_data_files.bl_corr_files:
        bl_corr_xml = ET.SubElement(mcobs_handler_init,"BLCorrelatorData")
        for filename in list(this_data_handler.raw_data_files.bl_corr_files)+list(this_data_handler.averaged_data_files.bl_corr_files)+list(this_data_handler.rotated_data_files.bl_corr_files):
            flist = filename.xml()
            bl_corr_xml.insert(1,flist)
    if this_data_handler.raw_data_files.bl_vev_files or this_data_handler.averaged_data_files.bl_vev_files:
        bl_vev_files_xml = ET.SubElement(mcobs_handler_init,"BLVEVData")
        #file list info
    if this_data_handler.raw_data_files.bin_files or this_data_handler.averaged_data_files.bin_files or this_data_handler.rotated_data_files.bin_files:
        bin_files_xml = ET.SubElement(mcobs_handler_init,"BinData")
        for filename in list(this_data_handler.raw_data_files.bin_files)+list(this_data_handler.averaged_data_files.bin_files)+list(this_data_handler.rotated_data_files.bin_files):
            ET.SubElement(bin_files_xml,"FileName").text = filename
    if this_data_handler.raw_data_files.sampling_files or this_data_handler.averaged_data_files.sampling_files or this_data_handler.rotated_data_files.sampling_files or additional_sampling_files:
        sampling_files_xml = ET.SubElement(mcobs_handler_init,"SamplingData")
        for filename in list(this_data_handler.raw_data_files.sampling_files)+list(this_data_handler.averaged_data_files.sampling_files)+list(this_data_handler.rotated_data_files.sampling_files)+additional_sampling_files:
            ET.SubElement(sampling_files_xml,"FileName").text = filename

    mcobs_init = sigmond.XMLHandler()
    mcobs_init.set_from_string(ET.tostring(mcobs_handler_init).decode())
    mcobs_get_handler = sigmond.MCObsGetHandler(mcobs_init, project_info.bins_info,project_info.sampling_info)
    mcobs_handler = sigmond.MCObsHandler(mcobs_get_handler, project_info.precompute)
    mcobs_handler.setSamplingBegin()

    return mcobs_handler, mcobs_get_handler

#take list of sigmond.MCEstimate and convert to csv
def estimates_to_csv( estimates, data_file):
    df = estimates_to_df( estimates )
    df.to_csv(data_file, index=False, header=True)

#take list of sigmond.MCEstimate and convert to pandas dataframe
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

#after rotation task is completed, takes the list of log files and retrieves the 
    #pivot quality information
def get_pivot_info(log_list):
    all_pivot_info = {}
    for file in log_list:
        if os.path.exists(file):
            rotation_log = sigmond_log.RotationLog(file)
            for i in range(rotation_log.num_rotations):
                pivot_info = {}
                pivot_info["metric null space check"] = rotation_log.metric_null_space_message(i)
                pivot_info["metric cond. init."] = rotation_log.metric_condition(i,False)
                pivot_info["metric cond. init."] = rotation_log.metric_condition(i)
                pivot_info["matrix cond. init."] = rotation_log.matrix_condition(i,False)
                pivot_info["matrix cond. init."] = rotation_log.matrix_condition(i)
                all_pivot_info[rotation_log.channel(i)] = pivot_info
    
    return all_pivot_info

#use the minimizer lmder in sigmond for fitting; very fast
def sigmond_fit( task_input, fitop, minimizer_configs, fit_configs, 
                mcobs_handler, sampling_info, subvev, logfile, delete_samplings = False, sh_priors={}, scat_info = []):

    this_fit_info = fit_info.FitInfo( 
        fitop,
        fit_configs['model'],
        fit_configs['tmin'],
        fit_configs['tmax'],
        subvev,
        fit_configs['ratio'],
        fit_configs['exclude_times'],
        fit_configs['noise_cutoff'],
        fit_configs['non_interacting_operators'],
        -1,
        -1,
        sim_fit = fit_configs['sim_fit'],
        initial_params = fit_configs['initial_params'],
        priors = fit_configs['priors'],
    )
    #minimizer info
    this_minimizer_info = sigmond_info.getMinimizerInfo(minimizer_configs)
    fit_options = {
        'minimizer_info': this_minimizer_info,
        'sampling_mode': sampling_info,
        'sh_priors': sh_priors,
        'scattering_fit_info': scat_info,
    }
    
    if not fit_configs['ratio'] and not fit_configs['sim_fit']: #and not vary
        fittype = "TemporalCorrelatorFit"
    elif fit_configs['ratio']:
        fittype = "TemporalCorrelatorInteractionRatioFit"
    elif fit_configs["sim_fit"]:
        fittype = "NSimTemporalCorrelatorFit/Fits/TemporalCorrelatorFit"

    task_input.doTemporalCorrelatorFit(this_fit_info,**fit_options) #add input for ni level priors

    task_xml = task_input.sigmondXML
    if delete_samplings:
        dummy_map = {}
        for i,(name,index) in enumerate(zip(task_xml.findall(f'TaskSequence/Task/{fittype}/Model/*/Name'),task_xml.findall(f'TaskSequence/Task/{fittype}/Model/*/IDIndex'))):
            if (name.text,index.text) not in dummy_map:
                dummy_map[(name.text,index.text)] = ("dummy", str(i))
                name.text = "dummy"
                index.text = str(i)
            else:
                name.text,  index.text = dummy_map[(name.text,index.text)]
            param = sigmond.MCObsInfo(name.text,int(index.text))
            mcobs_handler.eraseSamplings(param)
        if fit_configs['sim_fit']:
            task_xml.find(f'TaskSequence/Task/FinalEnergy/Name').text = "dummy"
            task_xml.find(f'TaskSequence/Task/FinalEnergy/IDIndex').text = "0"
            task_xml.find(f'TaskSequence/Task/FinalAmplitude/Name').text = "dummy"
            task_xml.find(f'TaskSequence/Task/FinalAmplitude/IDIndex').text = "1"
    
    if fit_configs['ratio']:
        #construct a "non-ratio" xml with ratio correlator as the fit correlator
        ratio_elem = task_xml.find(f'TaskSequence/Task/{fittype}/Ratio/*')
        optype = ratio_elem.tag
        opstring = ratio_elem.text

        fit_elem = task_xml.find(f'TaskSequence/Task/{fittype}')
        fit_elem.tag = "TemporalCorrelatorFit"
        ET.SubElement(fit_elem, optype).text = opstring

        #make the ratio correlator
        if optype == 'GIOperatorString':
            optype = sigmond.OpKind.GenIrrep
        else:
            optype = sigmond.OpKind.BasicLapH

        ni_ops = []
        for op in fit_configs['non_interacting_operators'].operators:
            ni_ops.append((op.operator_info, 0))
        ratio_op = sigmond.OperatorInfo(opstring, optype)
        sigmond.doCorrelatorInteractionRatioBySamplings(mcobs_handler,(this_fit_info.operator.operator_info, 0),ni_ops,
                                                        fit_configs['tmin'],fit_configs['tmax'],ratio_op, set(), 0)


    attempts = []
    if fit_configs['tmin_try_min'] and fit_configs['tmin_try_max'] and not delete_samplings:
        attempts+=list(range(fit_configs['tmin'],fit_configs['tmin_try_max']+1))
        attempts+=list(range(fit_configs['tmin']-1,fit_configs['tmin_try_min']-1,-1))
        
    if attempts:
        for attempt in attempts:
            if attempt!=fit_configs['tmin']:
                task_xml.find(f'TaskSequence/Task/{fittype}/MinimumTimeSeparation').text=str(attempt)
                old_obs_name = this_fit_info.obs_name
                this_fit_info.tmin = attempt
                obs_list = task_xml.findall(f'TaskSequence/Task/{fittype}/Model/*/Name')
                for obs in obs_list:
                    obs.text = obs.text.replace(old_obs_name, this_fit_info.obs_name)

            setuptaskhandler = sigmond.XMLHandler()
            setuptaskhandler.set_from_string(task_input.to_str())

            if fit_configs['sim_fit']:
                fit_xml_tag = "NSimTemporalCorrelatorFit"
            else:
                fit_xml_tag = "TemporalCorrelatorFit"
            fitxml = sigmond.XMLHandler(setuptaskhandler,fit_xml_tag)
            if fit_configs['sim_fit']:
                RTC = sigmond.NSimRealTemporalCorrelatorFit(fitxml,mcobs_handler,0)
            else:
                RTC = sigmond.RealTemporalCorrelatorFit(fitxml,mcobs_handler,0)
            chisqr = 0.0
            qual = 0.0
            log_xml = sigmond.XMLHandler("FitLogFile")

            try:
                best_fit_results = sigmond.doChiSquareFitting(RTC, this_minimizer_info, chisqr, qual, log_xml)
                if attempt!=fit_configs['tmin']:
                    # this_fit_info.tmin = attempt
                    logging.warning(f"Failed fit for {fitop} with tmin of {fit_configs['tmin']}. Using tmin={attempt} instead.")
                break
            except Exception as err:
                if not delete_samplings and attempt==attempts[-1]:
                    f = open(logfile,"w+")
                    f.write(log_xml.output(1))
                    f.close()
                    logging.warning(f"Failed fit for {fitop} with tmin of {fit_configs['tmin']}. Using tmin={attempt} instead.")
                    raise RuntimeError(err)
            
    else:
        setuptaskhandler = sigmond.XMLHandler()
        setuptaskhandler.set_from_string(task_input.to_str())
        if fit_configs['sim_fit']:
            fit_xml_tag = "NSimTemporalCorrelatorFit"
        else:
            fit_xml_tag = "TemporalCorrelatorFit"
        fitxml = sigmond.XMLHandler(setuptaskhandler,fit_xml_tag)
        if fit_configs['sim_fit']:
            RTC = sigmond.NSimRealTemporalCorrelatorFit(fitxml,mcobs_handler,0)
        else:
            RTC = sigmond.RealTemporalCorrelatorFit(fitxml,mcobs_handler,0)
        chisqr = 0.0
        qual = 0.0
        log_xml = sigmond.XMLHandler("FitLogFile")

        try:
            best_fit_results = sigmond.doChiSquareFitting(RTC, this_minimizer_info, chisqr, qual, log_xml)
        except Exception as err:
            if not delete_samplings:
                f = open(logfile,"w+")
                f.write(log_xml.output(1))
                f.close()
                logging.info(f"Fit Log written to '{logfile}'.")
            raise RuntimeError(err)
        
    log_xml.seek_unique("DegreesOfFreedom")
    dof = log_xml.get_text_content()
    log_xml.seek_unique("ChiSquarePerDof")
    chisqr = log_xml.get_text_content()
    log_xml.seek_unique("FitQuality")
    qual = log_xml.get_text_content()

    if not delete_samplings:
        f = open(logfile,"w+")
        f.write(log_xml.output(1))
        f.close()
        # logging.info(f"Logfile written to {logfile}.")

    return this_fit_info, best_fit_results, chisqr, qual, dof

#use python minimizers to fit; very slow, set up for parallel
def scipy_fit(fitop, minimizer_configs, fit_configs, 
                mcobs_handler, subvev, herm, Nt, 
                n_nodes, nsamplings, delete_samplings = False, sh_priors={}, scat_info = []):
    # start_time = time.time()
    this_fit_info = fit_info.FitInfo( 
        fitop,
        fit_configs['model'],
        fit_configs['tmin'],
        fit_configs['tmax'],
        subvev,
        fit_configs['ratio'],
        fit_configs['exclude_times'],
        fit_configs['noise_cutoff'],
        fit_configs['non_interacting_operators'],
        -1,
        -1,
        sim_fit = fit_configs['sim_fit'],
        initial_params = fit_configs['initial_params'],
        priors = fit_configs['priors'],
    )
    #minimizer info
    this_minimizer_info = sigmond_info.getMinimizerInfo(minimizer_configs)
    model = this_fit_info.model.sigmond_object(Nt)

    fit_op = fitop
    if this_fit_info.ratio:
        ni_ops = []
        for op in fit_configs['non_interacting_operators'].operators:
            ni_ops.append((op.operator_info, 0))
        fit_op = fitop.ratio_op
        sigmond.doCorrelatorInteractionRatioBySamplings(mcobs_handler,(this_fit_info.operator.operator_info, 0),ni_ops,
                                                        fit_configs['tmin'],fit_configs['tmax'],fit_op.operator_info, set(), 0)

    #if sim fit
    fit_ops = [fit_op]
    models = {fit_op: model}
    tranges = {}
    fit_infos = {fit_op: this_fit_info}
    parameter_map = {fit_op:list(range(this_fit_info.num_params))}
    if this_fit_info.sim_fit:
        for scat in scat_info:
            tranges[scat.operator] = list(range(scat.tmin, scat.tmax+1))
        if this_fit_info.model==fit_info.FitModel.DegTwoExpConspiracy:
            fit_ops.append(scat_info[0].operator)
            models[scat_info[0].operator] = scat_info[0].model.sigmond_object(Nt)
            fit_infos[scat_info[0].operator] = scat_info[0]
            fit_infos[scat_info[0].operator].priors = sh_priors[scat_info[0].obs_name]
            #SqrtGapToSecondEnergy==SqrtGapToSecondEnergy
            parameter_map[scat_info[0].operator] = [None]*scat_info[0].num_params
            parameter_map[scat_info[0].operator][scat_info[0].param_names.index("SqrtGapToSecondEnergy")] = this_fit_info.param_names.index("SqrtGapToSecondEnergy")
            index = this_fit_info.num_params
            for i in range(scat_info[0].num_params):
                if parameter_map[scat_info[0].operator][i]==None:
                    parameter_map[scat_info[0].operator][i] = index
                    index+=1
        if this_fit_info.model==fit_info.FitModel.TwoExpConspiracy:
            i0 = 0; i1 = 1
            shift1 = sh_priors[scat_info[i0].obs_name]["SqrtGapToSecondEnergy"]["Mean"]
            shift2 = np.sqrt(sh_priors[scat_info[i1].obs_name]["SqrtGapToSecondEnergy"]["Mean"]**2-shift1*shift1)
            if np.isnan(shift2):
                i0 = 1; i1 = 0
                shift1 = sh_priors[scat_info[i0].obs_name]["SqrtGapToSecondEnergy"]["Mean"]
                shift2 = np.sqrt(sh_priors[scat_info[i1].obs_name]["SqrtGapToSecondEnergy"]["Mean"]**2-shift1*shift1)

            fit_ops.append(scat_info[i0].operator)
            models[scat_info[i0].operator] = scat_info[i0].model.sigmond_object(Nt)
            fit_infos[scat_info[i0].operator] = scat_info[i0]
            fit_infos[scat_info[i0].operator].priors = sh_priors[scat_info[i0].obs_name]

            fit_ops.append(scat_info[i1].operator)
            fit_infos[scat_info[i1].operator] = scat_info[i1]
            fit_infos[scat_info[i1].operator].priors = sh_priors[scat_info[i1].obs_name]
            fit_infos[scat_info[i1].operator].model = fit_info.FitModel.TimeForwardTwoExponentialForCons
            fit_infos[scat_info[i1].operator].priors["SqrtGapToSecondEnergy"]["Mean"] = shift2
            fit_infos[scat_info[i1].operator].priors["SqrtGapToSecondEnergy"]["Error"] = 0.8*shift2
            models[scat_info[i1].operator] = fit_infos[scat_info[i1].operator].model.sigmond_object(Nt)

            #SqrtGapToSecondEnergy==SqrtGapToSecondEnergy
            parameter_map[scat_info[i0].operator] = [None]*scat_info[i0].num_params
            parameter_map[scat_info[i0].operator][scat_info[i0].param_names.index("SqrtGapToSecondEnergy")] = this_fit_info.param_names.index("SqrtGapToSecondEnergy")
            index = this_fit_info.num_params
            for i in range(scat_info[i0].num_params):
                if parameter_map[scat_info[i0].operator][i]==None:
                    parameter_map[scat_info[i0].operator][i] = index
                    index+=1
            
            #SqrtGapToSecondEnergy==SqrtGapToSecondEnergyShift
            #SqrtGapToThirdEnergy==SqrtGapToSecondEnergy
            parameter_map[scat_info[i1].operator] = [None]*scat_info[i1].num_params
            parameter_map[scat_info[i1].operator][scat_info[i1].param_names.index("SqrtGapToSecondEnergyShift")] = this_fit_info.param_names.index("SqrtGapToSecondEnergy")
            parameter_map[scat_info[i1].operator][scat_info[i1].param_names.index("SqrtGapToSecondEnergy")] = this_fit_info.param_names.index("SqrtGapToThirdEnergy")
            for i in range(scat_info[i1].num_params):
                if parameter_map[scat_info[i1].operator][i]==None:
                    parameter_map[scat_info[i1].operator][i] = index
                    index+=1


    if delete_samplings:
        for op in parameter_map:
            for i in parameter_map[op]:
                mcobs_handler.eraseSamplings(sigmond.MCObsInfo("dummy", i))

    attempts = []
    if fit_configs['tmin_try_min'] and fit_configs['tmin_try_max'] and not delete_samplings:
        attempts+=list(range(fit_configs['tmin'],fit_configs['tmin_try_max']+1))
        attempts+=list(range(fit_configs['tmin']-1,fit_configs['tmin_try_min']-1,-1))

    if attempts:
        for attempt in attempts:
            tranges[fit_op] = list(range(attempt, this_fit_info.tmax))
            res = complete_one_fit(mcobs_handler, tranges, fit_ops, models, this_minimizer_info, fit_infos, 
                                   parameter_map, n_nodes, nsamplings, herm, subvev)
    else:
        tranges[fit_op] = list(range(this_fit_info.tmin, this_fit_info.tmax+1))
        res = complete_one_fit(mcobs_handler, tranges, fit_ops, models, this_minimizer_info, fit_infos, 
                               parameter_map, n_nodes, nsamplings, herm, subvev, delete_samplings)

    dof = -len(res.x)
    for tvals in tranges.values():
        dof+=len(tvals)
    # dof = len(tvals)-len(res.x)
    chisqr = res.fun/dof
    pval = sigmond.getChiSquareFitQuality(dof, res.fun)

    bestfit_params=[]
    for ip in range(this_fit_info.num_params):
        bestfit_params.append(mcobs_handler.getEstimate(this_fit_info.fit_param_obs(ip)))

    if this_fit_info.sim_fit:
        fit_infos[scat_info[0].operator].priors = {}
        fit_infos[scat_info[0].operator].model = fit_info.FitModel.TimeForwardTwoExponential
        if this_fit_info.model==fit_info.FitModel.TwoExpConspiracy:
            fit_infos[scat_info[1].operator].priors = {}
            fit_infos[scat_info[1].operator].model = fit_info.FitModel.TimeForwardTwoExponential

    return this_fit_info, bestfit_params, chisqr, pval, dof

#scipy fits to one observable over all its samples
def complete_one_fit(mcobs_handler, tvals, fit_op, model, this_minimizer_info, this_fit_info, 
                     parameter_map, n_nodes, nsamplings, herm, subvev, delete_samplings=False):

    central_fit = fit_op[0]
    sh_fits = fit_op[1:]
    if delete_samplings:
        param_obs_keys = [sigmond.MCObsInfo("dummy",ip) for ip in range(this_fit_info[central_fit].num_params)]
    else:
        param_obs_keys = [this_fit_info[central_fit].fit_param_obs(ip) for ip in range(this_fit_info[central_fit].num_params)]
    for sh_fit in sh_fits:
        if delete_samplings:
            for index in parameter_map[sh_fit]:
                if sigmond.MCObsInfo("dummy",index) not in param_obs_keys:
                    param_obs_keys.append(sigmond.MCObsInfo("dummy",index))
        else:
            new_obs_name = this_fit_info[sh_fit].obs_name+"-"+this_fit_info[central_fit].obs_name
            for ip,index in enumerate(parameter_map[sh_fit]):
                if index==len(param_obs_keys):
                    param_obs_keys.append(sigmond.MCObsInfo(new_obs_name,this_fit_info[sh_fit].obs_id(ip)))
                elif index>len(param_obs_keys):
                    logging.error("Ive made a mistake")


    priors = {}        
    for prior in this_fit_info[central_fit].priors:
        if prior in this_fit_info[central_fit].param_names:
            priors[this_fit_info[central_fit].param_names.index(prior)] = this_fit_info[central_fit].priors[prior]

    for sh_fit in sh_fits:
        for prior in this_fit_info[sh_fit].priors:
            if prior in this_fit_info[sh_fit].param_names:
                new_prior_index = parameter_map[sh_fit][this_fit_info[sh_fit].param_names.index(prior)]
                print(sh_fit, prior, new_prior_index)
                if new_prior_index in priors:
                    logging.error("this is bad")
                priors[new_prior_index] = this_fit_info[sh_fit].priors[prior]

    #retrieve all samples and covmatrix
    ntvals = 0
    for op in tvals:
        ntvals += len(tvals[op])
    all_datapoints = np.zeros((nsamplings+1,ntvals))
    all_covmatrices = np.zeros((nsamplings+1,ntvals,ntvals))
    matrix_rows = []
    for op in fit_op:
        matrix_rows += [(op,t) for t in tvals[op]]
    mcobs_handler.setSamplingBegin()
    isamp = 0
    while not mcobs_handler.isSamplingEnd():
        old_op = None
        for i, (op, t) in enumerate(matrix_rows):
            if op!=old_op:
                corr = sigmond.CorrelatorAtTimeInfo(op.operator_info, op.operator_info, t, herm, subvev)
                old_op = op
            corr.resetTimeSeparation(t)
            obs = sigmond.MCObsInfo(corr, sigmond.ComplexArg.RealPart)
            all_datapoints[isamp][i] = mcobs_handler.getCurrentSamplingValue(obs)
            old_op2 = None
            for j, (op2, t2) in enumerate(matrix_rows[:i+1]):
                if op2!=old_op2:
                    corr2 = sigmond.CorrelatorAtTimeInfo(op2.operator_info, op2.operator_info, t2, herm, subvev)
                    old_op2 = op2
                corr2.resetTimeSeparation(t2)
                obs2 = sigmond.MCObsInfo(corr2, sigmond.ComplexArg.RealPart)
                # if op==op2:
                all_covmatrices[isamp,i,j] = mcobs_handler.getCovariance(obs, obs2)
                all_covmatrices[isamp,j,i] = np.conj(mcobs_handler.getCovariance(obs, obs2))

        mcobs_handler.setSamplingNext()
        isamp+=1

    mcobs_handler.setSamplingBegin()
    #fit mean
    res = minimize_sample(all_datapoints[0], all_covmatrices[0], tvals, fit_op, model, this_minimizer_info, 
                          parameter_map, herm, subvev, priors, global_fit=True)
    for ip, param in enumerate(res.x):
        mcobs_handler.putCurrentSamplingValue(param_obs_keys[ip], param, True)

    if not res.success:
        raise RuntimeError("Fit on the mean failed.")

    #fit samples
    threads = [None]*n_nodes
    thread_keys = [None]*n_nodes
    ithread = 0
    sampling_res = {}
    logging.info(f"\t\tDistributing samples on {n_nodes} threads...")
    for isamp in range(1,nsamplings+1):

        if not isamp%100:
            logging.info(f"\t\t{isamp} samplings computed on thread {ithread}...")

        sample_priors = {}
        sample_priors.update(priors)
        for ip in priors:
            sample_priors[ip]["Mean"] = np.random.normal(priors[ip]["Mean"],2.0*priors[ip]["Error"])

        if threads[ithread]!=None:
            threads[ithread].join()
            if not sampling_res[thread_keys[ithread]].success:
                for thread in threads:
                    if thread!=None:
                        thread.join()
                raise RuntimeError("Fit on the one of the samplings failed.")

        thread_keys[ithread] = isamp

        threads[ithread] = Thread(target=minimize_sample, args=(all_datapoints[isamp], all_covmatrices[isamp], tvals, fit_op, model, this_minimizer_info, 
                                parameter_map, herm, subvev, sample_priors, res.x, False, None, sampling_res, isamp))
        # sampling_res = minimize_sample(all_datapoints[isamp], all_covmatrices[isamp], tvals, fit_op, model, this_minimizer_info, 
        #                             parameter_map, herm, subvev, sample_priors, res.x)
        threads[ithread].start()
        ithread+=1
        ithread = ithread%n_nodes

    for thread in threads:
        if thread!=None:
            thread.join()

    #save samples
    mcobs_handler.setSamplingNext()
    isamp = 1
    while not mcobs_handler.isSamplingEnd():
        if not sampling_res[isamp].success:
            raise RuntimeError("Fit on the one of the samplings failed.")
        
        for ip, param in enumerate(sampling_res[isamp].x):
            mcobs_handler.putCurrentSamplingValue(param_obs_keys[ip], param, True)

        mcobs_handler.setSamplingNext()
        isamp+=1


    return res #fit results on the mean

#one scipy fit to one sample
def minimize_sample(datapoints, cov_matrix, tvals, fitop, model, this_minimizer_info, parameter_map, herm, subvev, 
                    priors, initial_params=[], global_fit=False, comm=None, results = {}, index = 0):

    if not len(initial_params):
        num_parameters = 0
        for op in fitop:
            num_parameters = max(num_parameters, max(parameter_map[op]))
        num_parameters += 1
        initial_params = [None]*num_parameters
        this_op_datapoints = None
        for i,op in enumerate(fitop):
            if i==0:
                this_op_datapoints = datapoints[:len(tvals[op])]
            else:
                this_op_datapoints = datapoints[len(tvals[fitop[i-1]]):len(tvals[fitop[i-1]])+len(tvals[op])]

            initial_params_set = model[op].guessInitialParamValuesPy(this_op_datapoints,tvals[op])
            for i, ip in enumerate(parameter_map[op]):
                if initial_params[ip]!=None:
                    initial_params[ip] += initial_params_set[i] 
                    initial_params[ip] /= 2.0
                else:
                    initial_params[ip] = initial_params_set[i] 
            for ip in priors:
                initial_params[ip] = priors[ip]["Mean"]

    if global_fit:
        initial_params[0] = 1.403727
        gres = sp.optimize.basinhopping(minimize_corr_function, initial_params, stepsize=0.1, target_accept_rate=0.1,
                                        seed=0, #niter_success = 10, 
                                 minimizer_kwargs={"method":"Nelder-Mead",
                                                   "args":(cov_matrix, model, tvals, datapoints, fitop, parameter_map, priors)})
        print(gres.x, gres.fun)
        initial_params = gres.x


        # print(np.linalg.eigvalsh(cov_matrix))
        # print(initial_params, minimize_corr_function(initial_params, cov_matrix, model, tvals, datapoints, fitop, parameter_map, priors))
    # return initial_params
    res = sp.optimize.minimize(minimize_corr_function, initial_params, #bounds=bounds,
                               method="Powell", #jac=jac_matrix, #'2-point', 
                               args=(cov_matrix, model, tvals, datapoints, fitop, parameter_map, priors), 
                               tol=this_minimizer_info.getChiSquareRelativeTolerance(),
                               options={"maxiter":this_minimizer_info.getMaximumIterations()})

    if comm!=None:
        comm.send(res)

    results[index] = res

    return res

# def jac_matrix(parameters, cov_matrix, model, trange, datapoints):
#     matrix = []
#     for param in parameters:
#         matrix.append([])
    
#     for t in trange:
#         gradt = model.evalGradientPy(parameters, t)
#         for ip, grad in enumerate(gradt):
#             matrix[ip] = grad

#     return matrix

#the function minimized in scipy fit
def minimize_corr_function(parameters, cov_matrix, model, trange, datapoints, fitop, parameter_map, priors={}):
    #get model points from parameters and model
    modelpoints = []
    divied_parameters = {}
    for op in parameter_map:
        divied_parameters[op] = []
        for ip in parameter_map[op]:
            divied_parameters[op].append(parameters[ip])

    for op in fitop:
        for t in trange[op]:
            modelpoints.append(model[op].eval(divied_parameters[op], t))

    prior_sum = 0.0
    for ip in priors:
        prior_sum += (parameters[ip] - priors[ip]["Mean"])*(parameters[ip] - priors[ip]["Mean"])/priors[ip]["Error"]/priors[ip]["Error"]

    # print(parameters, prior_sum, correlated_chisquare(datapoints, cov_matrix, np.array(modelpoints), prior_sum))
    return correlated_chisquare(datapoints, cov_matrix, np.array(modelpoints), prior_sum)

#the correlated chi-square of a model and datapoint set
def correlated_chisquare(data, cov_matrix, modelpoints, prior_sum):

    # cov = np.kron(data, data) #sp.stats.Covariance.from_diagonal(data)
    invcov = np.linalg.inv(cov_matrix)
    
    residuals = data-modelpoints
    covxres = np.sum(np.multiply(invcov, residuals), axis=1)
    chisquare = np.sum(residuals*covxres)+prior_sum

    return chisquare

#from a file, pivot_file, retrieve the pivot type of that file.
def get_pivot_type(pivot_file):
    ftype = sigmond.getFileID(pivot_file)
    if ftype==sigmond.FileType.SinglePivot_CN or ftype==sigmond.FileType.SinglePivot_RN:
        pivot_type = "single_pivot"
    elif ftype==sigmond.FileType.RollingPivot:
        pivot_type = "rolling_pivot"
    else:
        return None
    return sigmond_info.PivotType(pivot_type)

#create a sigmond.pivoter object for handling pivot operations and info
def setup_pivoter(pivot_type, pivot_file, channel, mcobs):
    #set up pivoter based on pivot file
    xmlinitiate_str = f'<Task><{pivot_type.name}Initiate><ReadPivotFromFile>'
    xmlinitiate_str += f'<PivotFileName>{pivot_file}[{repr(channel)}]</PivotFileName>'
    xmlinitiate_str += f'</ReadPivotFromFile></{pivot_type.name}Initiate></Task>'
    xmlinitiate = sigmond.XMLHandler()
    xmlinitiate.set_from_string(xmlinitiate_str)
    argsinitiate = sigmond.ArgsHandler(xmlinitiate)
    loghelper = sigmond.LogHelper()
    xmlout = sigmond.XMLHandler("LogFile")

    pivoter = sigmond.Pivot()
    pivoter.setType(pivot_type.name)
    pivoter.initiatePivotPython(mcobs, argsinitiate, loghelper)
    pivoter.checkInitiate(loghelper, xmlout)
    return pivoter

#based on the 'only' setting in task_configs, retrieve mom info from 
    #all items of the list and return in a list of ints
def get_selected_mom( task_configs):
    only_moms = []
    if 'only' in task_configs:
        for item in task_configs['only']:
            if item.startswith('PSQ='):
                only_moms.append(int(item.replace('PSQ=',"")))
            elif item.startswith('psq='):
                only_moms.append(int(item.replace('psq=',"")))
            else:
                only_moms.append(channel.Channel.CreateFromString(item).psq)
    only_moms = list(set(only_moms))
    return only_moms

#filer list based on only and omit settings, only overrides omit
def filter_channels( task_configs, channel_list):
    final_channels = []
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
        for channel in channel_list:
            if channel.psq in only_moms or str(channel) in only_channels:
                final_channels.append(channel)
            else:
                logging.info(f'Channel {str(channel)} omitted due to "only" setting.')
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
        for channel in channel_list:
            if channel.psq in omit_moms or str(channel) in omit_channels:
                logging.info(f'Channel {str(channel)} omitted due to "omit" setting.')
            else:
                final_channels.append(channel)
        # final_moms = list(set(range(10))-set(omit_moms))
    else: 
        final_channels = channel_list[:]
    return final_channels

#from estimates info, plot all correlators in operators matrix. If data is None, then
    #it assumed that the estimates are saved to file
def write_channel_plots(operators, plh, create_pickles, create_pdfs, pdh, data=None):
    saved_to_self = (data!=None)
    for op1 in operators:
        for op2 in operators:
            corr = sigmond.CorrelatorInfo(op1.operator_info,op2.operator_info)
            corr_name = repr(corr).replace(" ","-")

            plot_files = [
                pdh.corr_plot_file( corr_name, "pickle"),
                pdh.corr_plot_file( corr_name, "pdf"),
                pdh.effen_plot_file( corr_name, "pickle"),
                pdh.effen_plot_file(corr_name, "pdf"),
            ]
            for plot_file in plot_files:
                if os.path.isfile(plot_file):
                    os.remove(plot_file)

            try:
                if saved_to_self:
                    df = data[op1][op2]["corr"]
                else:
                    df = pd.read_csv(pdh.corr_estimates_file(corr_name))
            except pd.errors.EmptyDataError as err:
                logging.warning(f"pandas.errors.EmptyDataError: {err} for correlator {corr_name}.")

            plh.clf()
            try:
                plh.correlator_plot(df, 0) #, op1, op2) #0 for regular corr plot
            except KeyError as err:
                logging.warning(f"No correlator estimates could be optained for {corr_name}.")
                continue


            if create_pickles:
                plh.save_pickle(pdh.corr_plot_file( corr_name, "pickle"))
            if create_pdfs:
                plh.save_pdf(pdh.corr_plot_file( corr_name, "pdf"))

            try:
                if saved_to_self:
                    df = data[op1][op2]["effen"]
                else:
                    df = pd.read_csv(pdh.effen_estimates_file(corr_name))
                    
                plh.clf()
                plh.correlator_plot(df, 1) #, op1, op2) #1 for effective energy plot

                if create_pickles:
                    plh.save_pickle(pdh.effen_plot_file( corr_name, "pickle"))
                if create_pdfs:
                    plh.save_pdf( pdh.effen_plot_file(corr_name, "pdf")) 
            except pd.errors.EmptyDataError as err:
                # logging.warning(f"No effective energy estimates could be optained for {corr_name}.")
                pass
            except KeyError as err:
                # logging.warning(f"No effective energy estimates could be optained for {corr_name}.")
                pass
            
#sort channel based on isospin-strangeness-momentum
def channel_sort(item):
    return f"{item.isospin}{item.strangeness}{item.psq}"

#for processes, updates the iteration of the process index
    #such that the index ip is always less that nnodes
def update_process_index(ip,nnodes):
    ip += 1
    ip = ip%nnodes
    return ip