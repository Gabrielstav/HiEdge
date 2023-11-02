# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import yaml
from dataclasses import dataclass
from typing import List, Dict, NamedTuple, Optional
from pathlib import Path

class GenomicRange(NamedTuple):
    start: int
    end: int

@dataclass(frozen=True)
class ReferencePaths:
    blacklist_dir: Optional[Path] = None
    cytoband_dir: Optional[Path] = None

@dataclass(frozen=True)
class Paths:
    input_dir: Path
    output_dir: Path
    nchg_path: Path
    hg19: ReferencePaths
    hg38: ReferencePaths

@dataclass(frozen=True)
class PipelineSettings:
    reference_genome: str
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
    select_chromosomes: List[str]
    select_regions: Dict[str, List[GenomicRange]]
    omit_regions: Dict[str, List[GenomicRange]]
    fdr_threshold: float
    n_quantiles: int

    def __post_init__(self):
        object.__setattr__(self, 'select_regions', self._parse_ranges(self.select_regions))
        object.__setattr__(self, 'omit_regions', self._parse_ranges(self.omit_regions))

    @staticmethod
    def _parse_ranges(ranges_dict):
        parsed = {}
        for chrom, ranges in ranges_dict.items():
            range_list = []
            for rg in ranges.split(','):
                start, end = map(int, rg.split('-'))
                range_list.append(GenomicRange(start, end))
            parsed[chrom] = range_list
        return parsed

@dataclass(frozen=True)
class DaskSettings:
    chunksize: int
    work_stealing: bool
    default: str
    target: float
    spill: float
    pause: float
    terminate: float
    no_termination: bool


@dataclass(frozen=True)
class Config:
    version: float
    paths: Paths
    pipeline_settings: PipelineSettings
    dask_settings: DaskSettings


class ConfigMapper:
    def __init__(self, configuration_path):
        self.configuration_path = configuration_path
        self.config_data = self._load_configuration_from_file()

    def _load_configuration_from_file(self) -> Config:
        with open(self.configuration_path, "r") as config_file:
            config_data = yaml.safe_load(config_file)
            return Config(**config_data)



    
