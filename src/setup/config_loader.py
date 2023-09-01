# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import yaml

class ConfigLoader:
    def __init__(self, configuration_filepath):
        self.configuration_filepath = configuration_filepath
        self._loaded_config_data = None
        self._load_configuration_from_file()

    def _load_configuration_from_file(self):
        with open(self.configuration_filepath, 'r') as config_file:
            self._loaded_config_data = yaml.safe_load(config_file)

    def get_value_by_keys(self, *key_path):
        """Retrieve a configuration value using a key path."""
        current_data_level = self._loaded_config_data
        for key in key_path:
            if key in current_data_level:
                current_data_level = current_data_level[key]
            else:
                return None  # or raise an exception
        return current_data_level

def try_config_loader():
    config = ConfigLoader("/Users/GBS/PythonProjects/HiEdge/config/config.yaml")

    # Fetching a few values as an example
    input_dir = config.get_value_by_keys('paths', 'input_dir')
    output_dir = config.get_value_by_keys('paths', 'output_dir')
    threads = config.get_value_by_keys('settings', 'threads')
    fdr_threshold = config.get_value_by_keys('settings', 'fdr_threshold')

    print(f"Input Directory: {input_dir}")
    print(f"Output Directory: {output_dir}")
    print(f"Threads: {threads}")
    print(f"FDR Threshold: {fdr_threshold}")

    # Testing a non-existent key-path for graceful handling
    non_existent_value = config.get_value_by_keys('settings', 'non_existent_key')
    if non_existent_value is None:
        print("Successfully returned None for a non-existent key path!")

try_config_loader()



    
