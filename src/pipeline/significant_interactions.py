# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from src.setup.config_loader import Config
from src.setup.data_structures import SplineInput
import dask.dataframe as dd
from scipy.stats import binom


# TODO:
#   finish implementing these classes
#   need to handle no bias, inter and intra, normalization by all possible interactions and not normalizing
#   Just make classes for binomial, wallenius noncentral hypergeometric and fischer noncentral hypergeometric tests
#   that accepts the data and the splines (if intra --> nchg or binom) and only data (if inter --> fischer) and returns the p-values as new column in the dd.df
#   then make new classes for fdr and ci (for all data), these are called sequentially from a controller class with access to config, which is called from main pipeline
#   then make output format class 

class IntraPValueCalculator:

    def __init__(self, data, spline, total_interactions, config):
        self.data = data
        self.spline = spline
        self.total_interactions = total_interactions
        self.config = config

    def _calculate_expected_frequency(self):
        # Expected frequency for all intra interactions from the spline
        self.data["expected_frequency"] = self.data["genomic_distance"].map_partitions(
            lambda x: self.spline(x), meta=("x", "float64"))

    def _calculate_p_values_binom(self):
        # Calculate p-values using a binomial test
        self.data["p_value_binom"] = self.data.map_partitions(
            lambda df: binom.sf(
                df["interaction_count"] - 1,
                self.total_interactions,
                df["expected_frequency"] * df["bias_1"] * df["bias_2"]
            ),
            meta=("x", "float64")
        )

    def calculate_p_valus_nchg(self):
        # calculate p-values using the nchg test
        # use wallenius noncentral hypergeometric distribution:
        self.data["p_value_nchg"] = self.data.map_partitions(
            lambda df:


    def run(self):
        self._calculate_expected_frequency()
        self._calculate_p_values_binom()
        return self.data[["genomic_distance", "interaction_count", "expected_frequency", "p_value"]]


class FDRCalculator:

    def __init__(self, p_values):
        self.p_values = p_values

    def calculate_fdr(self):
        # use BH method, enable diff thresh for intra and inter?
        pass


class SignificanceCalculatorController:

    def __init__(self, config, intra_data, inter_data):
        self.config = config
        self.intra_data = intra_data
        self.inter_data = inter_data

    def run(self):
        # call the intra and inter p-value calculators
        # using the config, then calculate the FDR
        # then return data by  instantiating a dataclass with the results
        pass
