# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
import pybedtools as pbt

class ReferenceFileParser:

    def __init__(self, config: Config):
        self.config = config
        self.reference_genome = config.pipeline_settings.reference_genome
        self.blacklist_data = self.parse_blacklist_file()
        self.cytoband_data = self.parse_cytoband_file()

    def parse_blacklist_file(self):
        if self.config.pipeline_settings.reference_genome == "hg19":
            blacklist_path = self.config.paths.hg19.blacklist_dir
        elif self.config.pipeline_settings.reference_genome == "hg38":
            blacklist_path = self.config.paths.hg38.blacklist_dir
        else:
            raise ValueError(f"Unsupported reference genome: {self.reference_genome}")
        return pbt.BedTool(blacklist_path)

    def parse_cytoband_file(self):
        if self.config.pipeline_settings.reference_genome == "hg19":
            cytoband_path = self.config.paths.hg19.cytoband_dir
        elif self.config.pipeline_settings.reference_genome == "hg38":
            cytoband_path = self.config.paths.hg38.cytoband_dir
        else:
            raise ValueError(f"Unsupported rerefernec genome: {self.reference_genome}")
        return pbt.BedTool(cytoband_path)
