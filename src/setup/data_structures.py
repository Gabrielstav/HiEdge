# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from pathlib import Path
import dask.dataframe as dd

@dataclass
class Metadata:
    experiment: str
    resolution: int

@dataclass
class GroupedFiles:
    metadata: Metadata
    bed_file: Path
    matrix_file: Path

@dataclass
class BedpeOutput:
    metadata: Metadata
    bedpe_ddf: dd.DataFrame

@dataclass
class BlacklistOutput:
    metadata: Metadata
    blacklist_ddf: dd.DataFrame

@dataclass
class CytobandOutput:
    metadata: Metadata
    cytoband_ddf: dd.DataFrame

@dataclass
class NchgOutput:
    metadata: Metadata
    nchg_ddf: dd.DataFrame  # idk yet

@dataclass
class PadjOutput:
    metadata: Metadata
    padj_ddf: dd.DataFrame

@dataclass
class WriteEdgelist:
    metadata: Metadata
    edgelist: Path
