# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from src.setup.config_loader import Config
from src.setup.data_structures import SplineInput
import dask.dataframe as dd
from scipy.stats import binom
# use binom survival function instead? scsp.bdtrc


# TODO:
#   finish implementing these classes

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

    def _calculate_p_values(self):
        # Calculate p-values using a binomial test
        self.data["p_value"] = self.data.map_partitions(
            lambda df: binom.sf(
                df["interaction_count"] - 1,
                self.total_interactions,
                df["expected_frequency"] * df["bias_1"] * df["bias_2"]
            ),
            meta=("x", "float64")
        )

    def run(self):
        self._calculate_expected_frequency()
        self._calculate_p_values()
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
    