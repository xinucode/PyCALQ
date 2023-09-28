import general.config_handler as ch

class LuscherSchmuscher:

    def __init__( self, general_configs, task_configs = {} ):
        self.general_configs = ch.ConfigHandler(general_configs).configs
        self.task_configs = ch.ConfigHandler(task_configs).configs
        
        #reformat configs/import them from files, combine task dicts if multiple
        
        #perform checks on configs
        
        
    def run( self ):
        print("Hello World")
        print("general_configs:", self.general_configs)
        print("task_configs:", self.task_configs)
        #probably perform the tasks in an order that makes sense