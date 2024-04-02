import logging
import os
import xmltodict
import xml.etree.ElementTree as ET
import matplotlib
import shutil
from typing import NamedTuple
import pandas as pd

import sigmond
from sigmond_scripts.analysis.data_handling import data_files, data_handler
from sigmond_scripts.analysis.sigmond_info import fit_info, sigmond_info
from sigmond_scripts.analysis.operator_info import channel


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

    datafiles = data_files.DataFiles()

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
    mcobs_handler, mcobs_get_handler = get_mcobs_handlers(this_data_handler, project_info)
    return this_data_handler, mcobs_handler, mcobs_get_handler

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

    # filemap = sigmond.XMLHandler()
    # mcobs_get_handler.getFileMap(filemap)
    # print(filemap.output(0))

    mcobs_handler = sigmond.MCObsHandler(mcobs_get_handler, project_info.precompute)
    mcobs_handler.setSamplingBegin()

    return mcobs_handler, mcobs_get_handler


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
                # print(this_fit_info.model_name)
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
                    # print(old_obs_name, this_fit_info.obs_name, obs.text)

            setuptaskhandler = sigmond.XMLHandler()
            setuptaskhandler.set_from_string(task_input.to_str())

            if fit_configs['sim_fit']:
                # print(task_input.to_str())
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
        # print(task_input.to_str())
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

def get_pivot_type(pivot_file):
    ftype = sigmond.getFileID(pivot_file)
    if ftype==sigmond.FileType.SinglePivot_CN or ftype==sigmond.FileType.SinglePivot_RN:
        pivot_type = "single_pivot"
    elif ftype==sigmond.FileType.RollingPivot:
        pivot_type = "rolling_pivot"
    else:
        return None
    return sigmond_info.PivotType(pivot_type)

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
        final_moms = list(set(range(10))-set(omit_moms))
    else: 
        final_channels = channel_list[:]
    return final_channels


def write_channel_plots(operators, plh, create_pickles, create_pdfs, pdh, data=None):
    saved_to_self = (data!=None)
    for op1 in operators:
        for op2 in operators:
            corr = sigmond.CorrelatorInfo(op1.operator_info,op2.operator_info)
            corr_name = repr(corr).replace(" ","-")

            try:
                if saved_to_self:
                    df = data[op1][op2]["corr"]
                else:
                    df = pd.read_csv(pdh.corr_estimates_file(corr_name))
            except pd.errors.EmptyDataError as err:
                logging.warning(f"pandas.errors.EmptyDataError: {err} for correlator {corr_name}.")

            plh.clf()
            plh.correlator_plot(df, 0) #, op1, op2) #0 for regular corr plot

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
                pass
                # continue
                # logging.warning(f"pandas.errors.EmptyDataError: {err} for effective energy {corr_name}.")
            
#sort channel based on isospin-strangeness-momentum
def channel_sort(item):
    return f"{item.isospin}{item.strangeness}{item.psq}"
