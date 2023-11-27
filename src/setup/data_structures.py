# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from pathlib import Path
import dask.dataframe as dd
import dask.array as da
from typing import Optional

# rename file to pipeline_controleller from data_structures? and use nested dataclass where metadata is the same and fields are the structs for each pipeline method?
# could be nice to have metadata in one place since its not changing for one experiment, different instances of this class would have different experiment values
# so that each experiment is one instance, then withtin this experiemtn we have all the resolutions and data ?
# Could be potentially troublesome in regards to the config file, because the config dicitates the processing steps, but this change would more accurately reflect the
# actual structure of the data since each experiment can have many different resolutions, but we also want to be able to run many experiments in parallel with the same
# config if the user wants to (our file discovey and gropuing classes supports this functionality, we'd just need to change what field/class the metadata is instantiated to)?


@dataclass
class Metadata:
    experiment: str
    resolution: int
    interaction_type: Optional[str] = None

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
class SplittingOutput:
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


# IDK YET:
@dataclass
class NchgOutput:
    metadata: Metadata
    nchg_ddf: da.array

@dataclass
class PadjOutput:
    metadata: Metadata
    padj_ddf: dd.DataFrame

@dataclass
class WriteEdgelist:
    metadata: Metadata
    edgelist: Path
