import logging

import general.config_handler as ch
import fvspectrum.sigmond_view_corrs

# Thanks to Drew and https://stackoverflow.com/a/48201163/191474
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
            "view_corrs": None,
            # "average":{}
        }
}
TASK_ORDER = ["view_corrs"] #manage order of tasks
TASK_MAP = { #manage which classes to use for each unique task -> change for selection (fvspectrum)
    "view_corrs": fvspectrum.sigmond_view_corrs.SigmondViewCorrs
}


class LuscherSchmuscher:

    def __init__( self, general_configs, task_configs = DEFAULT_TASKS ):
        self.general_configs = ch.ConfigHandler(general_configs).configs['general']
        self.task_configs = ch.ConfigHandler(task_configs).configs['tasks']
        
        #perform checks on configs
            #require a project directory?

            #make sure that there is only one of each task
            #make sure that the tasks that are dependent on eachother have those
                #tasks -> discuss what to do in that case
        
    def run( self ):
        print("general_configs:", self.general_configs)
        print("task_configs:", self.task_configs)

        #output default input/default configurations somewhere in the project directory

        #probably perform the tasks in an order that makes sense
        for task in TASK_ORDER:
            if task in self.task_configs.keys():
                logging.info(f"Beginning task: {task}...")
                TASK_MAP[task](task, self.general_configs,self.task_configs[task])
                logging.info(f"Task {task} completed.")