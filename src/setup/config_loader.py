# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import yaml
from typing import List, Dict, Optional
from pathlib import Path
from pydantic import BaseModel, field_validator, Field


class GenomicRange(BaseModel):
    start: int
    end: int

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
    def log_dir(self):
        return self.run_dir / "logs"

    @property
    def tmp_dir(self):
        return self.run_dir / "tmp"


class PipelineSettings(BaseModel):
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
    remove_chromosomes: List[str]
    select_chromosomes: List[str]
    select_specific_regions: bool
    select_regions: Dict[str, List[str]] = Field(default_factory=dict)
    omit_regions: Dict[str, List[str]] = Field(default_factory=dict)
    use_interaction_distance_filters: bool
    interaction_distance_filters: Dict[int, InteractionDistanceFilter]
    output_format: str

    @field_validator("select_regions", "omit_regions")
    def convert_ranges_to_objects(cls, v: Dict[str, List[str]]):
        parsed = {}
        for chrom, ranges in v.items():
            range_objects = [GenomicRange(start=int(r.split("-")[0]), end=int(r.split("-")[1])) for r in ranges]
            parsed[chrom] = range_objects
        return parsed

    @field_validator("interactions_distance_filters", check_fields=False)
    def validate_interaction_distance_filters(cls, v: Dict[int, InteractionDistanceFilter]):
        parsed = {}
        for k, val in v.items():
            if not isinstance(k, int):
                raise ValueError(f"Keys in interaction_distance_filters must be integers, got {type(k).__name__}")
            if not isinstance(val, InteractionDistanceFilter):
                raise ValueError(f"Values in interaction_distance_filters must be InteractionDistanceFilter instances (integers), got {type(val).__name__}")
            parsed[k] = val


class StatisticalSettings(BaseModel):
    spline_passes: int
    fdr_threshold: float
    metabin_occupancy: int
    use_hicpro_bias: bool
    bias_lower_bound: float
    bias_upper_bound: float
    normalize_inter_expected_frequency: bool
    normalize_intra_expected_frequency: bool


class DaskSettings(BaseModel):
    chunksize: int
    work_stealing: bool
    default: str
    target: float
    spill: float
    pause: float
    terminate: float
    no_termination: bool


class Config(BaseModel):
    version: float
    run_name: str
    paths: Paths
    pipeline_settings: PipelineSettings
    statistical_settings: StatisticalSettings
    dask_settings: DaskSettings


class InstantiateConfig:
    def __init__(self, config_path: Path):
        self.config_path = config_path

    def load_configuration_from_file(self) -> Config:
        """Load and convert configuration from YAML file to Config class."""
        with open(self.config_path, "r") as config_file:
            config_data = yaml.safe_load(config_file)
            return Config(**config_data)
