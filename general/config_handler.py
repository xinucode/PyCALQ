import yaml
import json
import xmltodict
import os
import logging

class ConfigHandler:

    def __init__( self, input ):
        self.configs = {}
        self.init( input )

    def init( self, input ):
        if type(input)==str:
            yaml_file = False; json_file = False; xml_file = False
            try:
                with open(input, 'r') as f:
                    config_info = yaml.safe_load(f)
                    yaml_file = True
            except:
                logging.warning(f"{input} is not a YAML file")
                
            if not yaml_file:
                try:
                    with open(input, 'r') as f:
                        config_info = json.load(f)
                        json_file = True
                except:
                    logging.warning(f"{input} is not a JSON file")
                    
            if not yaml_file and not json_file:
                try:
                    with open(input, 'r') as f:
                        config_info = json.load(f)
                        xml_file = True
                except:
                    logging.warning(f"{input} is not a XML file")
                    
            if not yaml_file and not json_file and not xml_file:
                logging.error("Configs cannot be handled at this time.")
                
            self.configs.update(config_info)
                
        elif type(input)==dict:
            self.configs.update(input)
            
        elif type(input)==list:
            for this_input in input:
                self.init(this_input )
        else:
            logging.error("Configs cannot be handled at this time.")
