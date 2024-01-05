# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from pathlib import Path
from dataclasses import dataclass
import dask.dataframe as dd
from src.setup.config_loader import Config

# TODO:
#   This isdone after spline fitting, when calculating P values.
#   since bias is locus specific, we need to read in the bias into a dask df or other DS (array, dict) and then look up the bias for each locus
#   if a locus has a bias below or above a certain threshold, we need to filter it out, meaning all interactions that contains that locus will be removed.
#   So we need to:
#       1. read in the bias file
#       2. parse it into a DS
#       3. look up the bias for each locus
#       4. filter out interactions that contain loci with bias above or below a certain threshold
#       5. return the filtered interactions

class BiasApplier:

    # Apply bias correction using the bias_reader class (not yet implemented).

    def __init__(self, config: Config, bias_data=None):
        self.config = config
        self.bias_data = bias_data

    @staticmethod
    def apply_bias(data: dd.DataFrame) -> dd.DataFrame:
        # Apply bias correction here
        # Make sure to handle NaN and invalid values appropriately
        return data.persist()