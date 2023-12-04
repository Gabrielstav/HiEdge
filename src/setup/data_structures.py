# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from pathlib import Path
import dask.dataframe as dd
import dask.array as da
from typing import Optional

@dataclass
class Metadata:
    experiment: str
    resolution: int
    interaction_type: Optional[str] = None
    bias_file_path: Optional[Path] = None

@dataclass
class GroupedFiles:
    metadata: Metadata
    bed_file: Path
    matrix_file: Path

@dataclass
class BedpeOutput:
    data: dd.DataFrame
    metadata: Metadata

@dataclass
class FilteringOutput:
    data: dd.DataFrame
    metadata: Metadata

@dataclass
class BlacklistOutput:
    metadata: Metadata
    bedpe_ddf: dd.DataFrame

@dataclass
class CytobandOutput:
    metadata: Metadata
    bedpe_ddf: dd.DataFrame

@dataclass
class StatInput:
    metadata: Metadata
    bedpe_ddf: dd.DataFrame
