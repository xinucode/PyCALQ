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
            try: #tested
                with open(input, 'r') as f:
                    config_info = yaml.safe_load(f)
                    yaml_file = True
            except:
                logging.warning(f"{input} is not a YAML file")
                
            if not yaml_file: #tested
                try:
                    with open(input, 'r') as f:
                        config_info = json.load(f)
                        json_file = True
                except:
                    logging.warning(f"{input} is not a JSON file")
                    
            if not yaml_file and not json_file: #tested
                try:
                    with open(input, 'r') as f:
                        config_info = xmltodict.parse(f.read())
                        xml_file = True
                        
                    #first level guaranteed to only have one index
                    config_info[list(config_info.keys())[0]] = dict(config_info[list(config_info.keys())[0]])
                    
                    #issue with xml everything that is input is interpreted as a string. need to do more input processing if taken seriously
                    
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

    def check(self, root, requirements ):

        #check that root is singular and that is the set one
        if len(self.configs.keys())>1:
            logging.error(f"Config cannot have more than one root. Only '{root}' is allowed.")

        if list(self.configs.keys())[0]!=root:
            logging.error(f"{root} config must have '{root}' as the root.")

        #check that required input keys are present
        for requirement in requirements:
            if type(requirement)!=dict:
                if requirement not in self.configs[root].keys():
                    logging.error(f"Missing parameter '{requirement}' in {root} config.")
            else:
                requirement_key = list(requirement.keys())[0]
                if requirement_key not in self.configs[root].keys():
                    logging.error(f"Missing parameter '{requirement_key}' in {root} config.")
                for requirementee in requirement[requirement_key]:
                    if requirementee not in self.configs[root][requirement_key].keys():
                        logging.error(f"Missing parameter '{requirementee}' in {root} config.")


