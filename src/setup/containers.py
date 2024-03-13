# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from pathlib import Path
import dask.dataframe as dd
from typing import Optional

@dataclass
class Metadata:
    experiment: str
    resolution: int
    interaction_type: Optional[str] = None
    bias_file_path: Optional[Path] = None
    max_possible_interacting_bins_intra: Optional[int] = None
    max_possible_interacting_bins_inter: Optional[int] = None
    max_possible_interacting_bins_per_chromosome_intra: Optional[dict] = None
    unique_bins_per_chromosome: Optional[dict] = None

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
class StatisticalOutput:
    metadata: Metadata
    data: dd.DataFrame

