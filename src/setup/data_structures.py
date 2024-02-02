# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from pathlib import Path
import dask.dataframe as dd
from typing import Optional

# TODO:
#   use schame to validate dd dataframes of interactions later?
@dataclass
class Interaction:
    chr1: str
    start1: int
    end1: int
    chr2: str
    start2: int
    end2: int
    midpoint_1: int
    midpoint_2: int
    bias_1: Optional[float] = None
    bias_2: Optional[float] = None
    genomic_distance: Optional[int] = None
    expected_frequency: Optional[float] = None
    p_value: Optional[float] = None

# Dataset metadata (cell line)
@dataclass
class Metadata:
    experiment: str
    resolution: int
    interaction_type: Optional[str] = None
    bias_file_path: Optional[Path] = None
    max_possible_interaction_count_intra: Optional[int] = None
    max_possible_interaction_count_inter: Optional[int] = None
    interaction_count_per_chromosome: Optional[dict] = None


# TODO:
#   Maybe just reinstantiate the InteractionOutput class with data instead of having different dataclass containers that are at the same abstraction lvl?:
@dataclass
class GroupedFiles:
    metadata: Metadata
    bed_file: Path
    matrix_file: Path

@dataclass
class InteractionOutput:
    metadata: Metadata
    data: dd.DataFrame

@dataclass
class FilteringOutput:
    metadata: Metadata
    data: dd.DataFrame

@dataclass
class SplineInput:
    metadata: Metadata
    data: dd.DataFrame

@dataclass
class StatisticsOutput:
    metadata: Metadata
    data: dd.DataFrame

