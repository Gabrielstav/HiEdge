# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import yaml
from dataclasses import dataclass
from typing import List
from pathlib import Path


@dataclass(frozen=True)
class Paths:
    input_dir: Path
    output_dir: Path
    nchg_path: Path
    blacklist_dir: Path
    cytoband_dir: Path

@dataclass(frozen=True)
class Settings:
    hicpro_raw_dirname: str
    hicpro_norm_dirname: str
    inter: bool
    intra: bool
    mixed: bool
    threads: int
    executor: str
    no_split: bool
    normalized_data: bool
    resolutions: List[int]
    filter_blacklist: bool
    filter_cytobands: bool
    fdr_threshold: float
    n_quantiles: int

@dataclass(frozen=True)
class Config:
    version: float
    paths: Paths
    settings: Settings


class ConfigMapper:
    def __init__(self, configuration_path):
        self.configuration_path = configuration_path
        self.config_data = self._load_configuration_from_file()

    def _load_configuration_from_file(self) -> Config:
        with open(self.configuration_path, "r") as config_file:
            config_data = yaml.safe_load(config_file)
            return Config(**config_data)



    
