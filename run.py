import argparse
import pycalq

#pick up the general and any task config files passed in arguments
def pickup_configs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', "--general", help="general configuration file")
    parser.add_argument('-t', "--tasks", nargs='+', help="task(s) configuration file(s)", required=False)
    args = parser.parse_args()
    config_file = args.general
    task_configs = args.tasks
    return config_file, task_configs

#runs pycalq
if __name__=="__main__":
    gen_configs, task_configs = pickup_configs()
    if task_configs:
        this_channel = pycalq.PyCALQ(gen_configs, task_configs)
    else:
        this_channel = pycalq.PyCALQ(gen_configs)
    this_channel.run()