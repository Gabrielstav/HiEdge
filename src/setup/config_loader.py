# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import yaml
from typing import List, Dict, Optional
from pathlib import Path
from pydantic import BaseModel, field_validator
from pydantic.dataclasses import dataclass


class GenomicRange(BaseModel):
    start: int
    end: int

    @staticmethod
    def from_string(s: str):
        start, end = map(int, s.split("-"))
        return GenomicRange(start=start, end=end)

class InteractionDistanceFilter(BaseModel):
    min_distance: int
    max_distance: int

class Version(BaseModel):
    version: float

class RunName(BaseModel):
    run_name: str

class ReferencePaths(BaseModel):
    blacklist_dir: Optional[Path] = None
    cytoband_dir: Optional[Path] = None

class Paths(BaseModel):
    input_dir: Path
    run_dir: Path
    hg19: ReferencePaths
    hg38: ReferencePaths

    @property
    def output_dir(self):
        return self.run_dir / "output"

    @property
    def plot_dir(self):
        return self.run_dir / "plots"

    @property
    def log_dir(self):
        return self.run_dir / "logs"

    @property
    def tmp_dir(self):
        return self.run_dir / "tmp"

@dataclass
class PipelineSettings:
    reference_genome: str
    hicpro_raw_dirname: str
    hicpro_norm_dirname: str
    interaction_type: str
    iced_data: bool
    round_iced_matrices: bool
    intra_resolutions: List[int]
    inter_resolutions: List[int]
    filter_blacklist: bool
    filter_cytobands: bool
    filter_self_interactions: bool
    remove_chromosomes: List[str]
    select_chromosomes: List[str]
    make_plots: bool
    select_specific_regions: bool
    omit_specific_regions: bool
    use_interaction_distance_filters: bool
    interaction_distance_filters: Dict[int, InteractionDistanceFilter]
    output_format: str
    output_type : str
    select_regions: Optional[Dict[str, List[str]]] = None
    omit_regions: Optional[Dict[str, List[str]]] = None

    def __post_init__(self):
        if self.select_regions is not None:
            self.select_regions = self.convert_ranges(self.select_regions)
        if self.omit_regions is not None:
            self.omit_regions = self.convert_ranges(self.omit_regions)

    @staticmethod
    def convert_ranges(regions: Dict[str, List[str]]) -> Dict[str, List[GenomicRange]]:
        return {
            chrom: [GenomicRange.from_string(range_str) for range_str in ranges]
            for chrom, ranges in regions.items()
        }

    @field_validator("interactions_distance_filters", check_fields=False)
    def validate_interaction_distance_filters(cls, v: Dict[int, InteractionDistanceFilter]):
        parsed = {}
        for k, val in v.items():
            if not isinstance(k, int):
                raise ValueError(f"Keys in interaction_distance_filters must be integers, got {type(k).__name__}")
            parsed[k] = val

class StatisticalSettings(BaseModel):
    spline_passes: int
    fdr_threshold: float
    metabin_occupancy: int
    use_hicpro_bias: bool
    bias_lower_bound: float
    bias_upper_bound: float
    use_filtered_data_for_average_contact_probability: bool
    use_sequential_fdr: bool


class Config(BaseModel):
    version: float
    run_name: str
    paths: Paths
    pipeline_settings: PipelineSettings
    statistical_settings: StatisticalSettings

class InstantiateConfig:
    def __init__(self, config_path: Path):
        self.config_path = config_path

    def load_configuration_from_file(self) -> Config:
        """Load and convert configuration from YAML file to Config class."""
        with open(self.config_path, "r") as config_file:
            config_data = yaml.safe_load(config_file)
            return Config(**config_data)
