# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import ConfigMapper as Config
from pathlib import Path

class BlacklistFileParser:

    def __init__(self, config: Config):
        self.config = config

    def parse_blacklist_file(self, blacklist_path: Path):


    def set_reference_genome(self):

        reference_genome = ""

        if self.config.config_data.pipeline_settings.reference_genome == "hg19":
            reference_genome = "hg19"
        if self.config.config_data.pipeline_settings.reference_genome == "hg38":
            reference_genome = "hg38"
        else:
            print(TypeError("Valid reference genome not found, needs to be hg19 or hg38"))

        return reference_genome
