# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import yaml
from dataclasses import dataclass, field
from typing import List, Dict, NamedTuple, Optional, get_origin, Type, Any, get_args
from pathlib import Path
import dataclasses

class GenomicRange(NamedTuple):
    start: int
    end: int

@dataclass
class Version:
    version: float

@dataclass
class RunName:
    run_name: str

@dataclass(frozen=True)
class ReferencePaths:
    blacklist_dir: Optional[Path] = None
    cytoband_dir: Optional[Path] = None

@dataclass(frozen=True)
class Paths:
    input_dir: Path
    run_dir: Path
    output_dir: Path = field(init=False)
    log_dir: Path = field(init=False)
    tmp_dir: Path = field(init=False)
    hg19: ReferencePaths
    hg38: ReferencePaths

    def __post_init__(self):
        object.__setattr__(self, "output_dir", self.run_dir / "output")
        object.__setattr__(self, "log_dir", self.run_dir / "logs")
        object.__setattr__(self, "tmp_dir", self.run_dir / "tmp")

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
    iced_data: bool
    max_interaction_range: int
    min_interaction_range: int
    resolutions: List[int]
    filter_blacklist: bool
    filter_cytobands: bool
    remove_chromosomes: List[str]
    select_chromosomes: List[str]
    select_regions: Dict[str, List[GenomicRange]]
    omit_regions: Dict[str, List[GenomicRange]]
    fdr_threshold: float
    metabin_occupancy: int
    spline_passes: int
    bias_lower_bound: float
    bias_upper_bound: float
    use_hicpro_bias: bool
    interaction_distance_filters: Dict[int, Dict[str, int]]
    normalize_inter_expected_frequency: bool
    normalize_intra_expected_frequency: bool
    output_format: str

    def __post_init__(self):
        object.__setattr__(self, "select_regions", self._parse_ranges(self.select_regions))
        object.__setattr__(self, "omit_regions", self._parse_ranges(self.omit_regions))
        object.__setattr__(self, "interaction_distance_filters", self._parse_distance_filters(self.interaction_distance_filters))

    @staticmethod
    def _parse_ranges(ranges_dict):
        parsed = {}
        for chrom, ranges in ranges_dict.items():
            range_list = []
            for rg in ranges.split(","):
                start, end = map(int, rg.split("-"))
                range_list.append(GenomicRange(start, end))
            parsed[chrom] = range_list
        return parsed

    @staticmethod
    def _parse_distance_filters(filters_dict):
        parsed = {}
        for resolution, ranges in filters_dict.items():
            parsed_resolution = {key: int(value) for key, value in ranges.items()}
            parsed[resolution] = parsed_resolution
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
    run_name: str
    paths: Paths
    pipeline_settings: PipelineSettings
    dask_settings: DaskSettings

def parse_genomic_ranges(range_str: str) -> List[GenomicRange]:
    """Parse a string of genomic ranges into a list of GenomicRange tuples."""
    return [GenomicRange(*map(int, part.split("-"))) for part in range_str.split(",")]

def convert_to_dataclass(target_type: Type[Any], data: Any) -> Any:
    """Convert dictionaries to dataclass instances, handling special types and lists."""
    if dataclasses.is_dataclass(target_type):
        field_values = {}
        for fields in dataclasses.fields(target_type):
            if not fields.init:  # Skip fields that are not initialized through the constructor
                continue
            if fields.name in data:
                field_data = data[fields.name]
                field_values[fields.name] = convert_field_value(fields.type, field_data)
            elif fields.default is not dataclasses.MISSING:
                field_values[fields.name] = fields.default
            elif fields.default_factory is not dataclasses.MISSING:
                field_values[fields.name] = fields.default_factory()
            else:
                continue  # Skip raising an error for missing fields with init=False
        return target_type(**field_values)
    else:
        return data

def convert_field_value(field_type: Type[Any], value: Any) -> Any:
    """Handle conversion for specific field types, including lists and dataclasses."""

    if dataclasses.is_dataclass(field_type):
        return convert_to_dataclass(field_type, value)

    elif field_type == Path:
        return Path(value)

    elif get_origin(field_type) == list:
        element_type = get_args(field_type)[0]
        if isinstance(value, str):
            # Handle special case for strings that need to be split and converted
            if element_type == GenomicRange:
                return parse_genomic_ranges(value)
            else:
                return [convert_field_value(element_type, v) for v in value.split(",")]
        else:
            return [convert_field_value(element_type, v) for v in value]
    else:
        # direct conversion for simple types or if no conversion is necessary
        return value

class ConfigMapper:
    def __init__(self, configuration_path: Path):
        self.configuration_path = configuration_path

    def load_configuration_from_file(self) -> Config:
        """Load and convert configuration from YAML file to Config dataclass."""
        with open(self.configuration_path, "r") as config_file:
            config_data = yaml.safe_load(config_file)
            return convert_to_dataclass(Config, config_data)
