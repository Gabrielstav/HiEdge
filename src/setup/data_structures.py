# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from pathlib import Path
import dask.dataframe as dd
import dask.array as da

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
    nchg_ddf: da.array  # ?? idk yet, probably dask.array since we need to vectorize data but no, we do not need to pass that to the downstream methods

@dataclass
class PadjOutput:
    metadata: Metadata
    padj_ddf: dd.DataFrame

@dataclass
class WriteEdgelist:
    metadata: Metadata
    edgelist: Path
