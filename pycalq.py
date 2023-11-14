import logging

import general.config_handler as ch
import general.project_directory as pd
import fvspectrum.sigmond_view_corrs

# Thanks to Drew and https://stackoverflow.com/a/48201163/191474
#ends code when run logging.error(message) or logging.critical(message)
# logging.warning(message) and logging.debug(message) won't end code but will output at certain verbosity settings
class ExitOnExceptionHandler(logging.StreamHandler):

  def emit(self, record):
    super().emit(record)
    if record.levelno in (logging.ERROR, logging.CRITICAL):
      raise SystemExit(-1)
    
logging.basicConfig(format='%(levelname)s: %(message)s', handlers=[ExitOnExceptionHandler()], level=logging.INFO)

# set up the tasks config and order and classes
DEFAULT_TASKS = { #manage default configurations
    "tasks":
        {
            "preview_corrs": None,
            # "average":{}
        }
}
                    
#manage order of tasks                                 
TASK_ORDER = ["preview_corrs", "average_corrs","rotate_corrs","fit_corrs",#correlator analysis  
              # leuscher analysis
              # includes load data, can plot mean data for check "data_load"
              # includes identifying channel thresholds in relevant energy region
              # singlew_channel_fit_mean seems good, need to idenfify single channels
              # signe channel fit mean can use bootstrap to estimate errors
              "single_channel_fit_mean"] #default is estimate errors, bootstrap option
              # "single_channel_fit_err?","coupled_channel_fit"] #luscher qc
TASK_MAP = { #manage which classes to use for each unique task -> change for selection (fvspectrum)
    "preview_corrs": fvspectrum.sigmond_view_corrs.SigmondPreviewCorrs,
}
TASK_DOC = { #imports documentation from each task
    "preview_corrs": fvspectrum.sigmond_view_corrs.doc,
}

#set required general parameters 
#items in list must be str or {str: list of str}
#we don't really want to go deeper than two levels
REQUIRED_GENERAL_CONFIGS = [
   'project_dir',
   'ensemble_id',
#    {'ensemble_info':['ensemble_id']},
]


class PyCALQ:

    def __init__( self, general_configs, task_configs = DEFAULT_TASKS ):
        self.general_configs_handler = ch.ConfigHandler(general_configs)
        self.task_configs_handler = ch.ConfigHandler(task_configs)

        #perform checks on configs
        self.general_configs_handler.check('general',REQUIRED_GENERAL_CONFIGS)
        self.task_configs_handler.check('tasks',[])
        
        self.general_configs = self.general_configs_handler.configs['general']
        self.task_configs = self.task_configs_handler.configs['tasks']

        self.proj_dir = pd.ProjectDirectoryHandler(self.general_configs['project_dir'])

    def run( self ):
        #output default input/default configurations somewhere in the project directory

        #probably perform the tasks in an order that makes sense
        for task in TASK_ORDER:
            if task in self.task_configs.keys():
                logging.info(f"Setting up task: {task}...")
                self.proj_dir.set_task(TASK_ORDER.index(task), task)
                this_task = TASK_MAP[task](task, self.proj_dir, self.general_configs,self.task_configs[task]) #initialize
                logging.info(f"Task {task} set up.")

                if 'info' in self.task_configs[task]:
                   if self.task_configs[task]['info']:
                      print(this_task.info)

                logging.info(f"Running task: {task}...")
                this_task.run() #perform the task, produce the data
                logging.info(f"Task {task} completed.")

                #do not plot if "plot: False" is present in task yaml
                if 'plot' in self.task_configs[task]:
                   if not self.task_configs[task]['plot']:
                      break
                   
                logging.info(f"Plotting task: {task}...")
                this_task.plot() 
                logging.info(f"Task {task} plotted.")




