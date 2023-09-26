import argparse
import luscher

# this is the script that runs everything else

# this script will ideally run one or some of the subset of tasks:
#  - compute the fv spectrum (Sarah, Arjun)
#  - use the luscher formalism to compute observables (Joseph)

# we will hopefully break this down and detail out these steps further in the future


def pickup_configs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', "--general", help="general configuration file")
    parser.add_argument('-t', "--tasks", nargs='+', help="task(s) configuration file(s)", required=False)
    args = parser.parse_args()
    config_file = args.general
    task_configs = args.tasks
    return config_file, task_configs

        
if __name__=="__main__":
    gen_configs, task_configs = pickup_configs()
    this_channel = luscher.LuscherSchmuscher(gen_configs, task_configs)
    this_channel.run()