# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from pathlib import Path
from dataclasses import dataclass
import dask.dataframe as dd
from src.setup.config_loader import Config

# TODO:
#   So we need to make two classes I think, one for reading in the bias, and one for applying it:
#       First class that reads in the bias file, can be stored here:
#       1. read in the bias file for each dataset (path to correct file stored in metadata field in dataclass containing data)
#       2. parse it into a dask df
#       3. calculate midpoints for each bias bin, set out-of-bounds bias to -1 (bounds are read from config instance)
#       Next class or method in the same class(?):
#       4. look up the bias for each locus, and sum the bias for each interaction. If either bias is -1, we set the interaction bias to -1.
#       4. filter out interactions that contain loci with bias -1, that is, remove them from the dataset.
#       5. return the filtered dataset, with the interaction bias added as a column.


class BiasReader:

    def __init__(self, config: Config, bias_file_path):
        self.config = config
        self.bias_file_path = bias_file_path

    def run_bias_reader(self):
        bias_df = self._read_bias()
        self._bias_quantiles(bias_df)  # log this later
        self._calculate_midpoints(bias_df)
        return bias_df

    def _read_bias(self):
        # Read bias file into Dask DataFrame
        bias_df = dd.read_csv(self.bias_file_path)
        return bias_df

    def _bias_quantiles(self, bias_df):
        quantiles = bias_df["bias"].quantile([0.05, 0.5, 0.95]).compute()
        # Log these in loggin output at some later stage
        return quantiles

    def _calculate_midpoints(self, bias_df):
        bias_df["midpoint"] = (bias_df["start"] + bias_df["end"]) / 2
        # Set out-of-bounds bias to -1
        bias_df.loc[bias_df["midpoint"] < self.config.pipeline_settings.bias_lower_bound, "bias"] = -1
        bias_df.loc[bias_df["midpoint"] > self.config.pipeline_settings.bias_upper_bound, "bias"] = -1
        bias_df.loc[bias_df["midpoint"] == "nan", "bias"] = -1
        return bias_df

class BiasApplier:

    def __init__(self, bias_data, config: Config):
        self.config = config
        self.bias_data = bias_data

    def apply_bias(self):


