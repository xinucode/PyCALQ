import argparse
# this is the script that runs everything else

# this script will ideally run one or some of the subset of tasks:
#  - compute the fv spectrum (Sarah, Arjun)
#  - use the luscher formalism to compute observables (Joseph)

# we will hopefully break this down and detail out these steps further in the future


def pickup_configs():
    parser = argparse.ArgumentParser()
    parser.add_argument("general_configs", help="general configuration file")
    parser.add_argument("task_configs", nargs='+', help="task(s) configuration file(s)")
    args = parser.parse_args()
    config_file = args.general_configs
    task_configs = args.task_configs
    return config_file, task_configs

class LuscherSchmuscher:

    def __init__( self, general_configs, task_configs ):
        self.general_configs = general_configs
        self.task_configs = task_configs
        
        #reformat configs/import them from files, combine task dicts if multiple
        
        #perform checks on configs
        
        
    def run( self ):
        print("Hello World")
        print("general_configs:", self.general_configs)
        print("task_configs:", self.task_configs)
        #probably perform the tasks in an order that makes sense
        
        
if __name__=="__main__":
    gen_configs, task_configs = pickup_configs()
    this_channel = LuscherSchmuscher(gen_configs, task_configs)
    this_channel.run()