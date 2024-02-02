# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from dataclasses import dataclass
from src.setup.config_loader import Config
from src.setup.data_structures import SplineInput
import dask.dataframe as dd
from scipy.stats import binom
from scipy.stats import nchypergeom_wallenius as nchg
from scipy.stats import binom


# TODO:
#  later: implement nchg test for intra p-vals

class IntraPValueCalculator:
    def __init__(self, data, spline, total_interactions, config: Config):
        self.data = data
        self.spline = spline
        self.total_interactions = total_interactions
        self.config = config

    def run(self):
        self._calculate_expected_frequency()
        self._calculate_p_values()
        return self.data[["genomic_distance", "interaction_count", "expected_frequency", "p_value"]]

    def _calculate_expected_frequency(self):
        # Calculate expected frequency using the spline
        self.data["expected_frequency"] = self.data["genomic_distance"].map_partitions(
            lambda x: self.spline(x), meta=("x", "float64"))

        if self.config.pipeline_settings.normalize_intra_expected_frequency:
            self.data["expected_frequency"] /= self.total_interactions

    def _calculate_p_values(self):
        # Check if bias should be used
        if self.config.pipeline_settings.use_hicpro_bias:
            bias_term = self.data["bias_1"] * self.data["bias_2"]
        else:
            bias_term = 1

        # Calculate p-values using a binomial test
        self.data["p_value"] = self.data.map_partitions(
            lambda df: binom.sf(
                df["interaction_count"] - 1,
                self.total_interactions,
                df["expected_frequency"] * bias_term
            ),
            meta=("x", "float64")
        )

class InterPValueCalculator:
    def __init__(self, config: Config, data, metadata, inter_prob, total_possible_interactions):
        self.config = config
        self.data = data
        self.metadata = metadata
        self.interChrProb = inter_prob
        self.total_possible_interactions = total_possible_interactions

    def run(self):
        self._calculate_prior_probability()
        self._calculate_p_values()
        return self.data[["interaction_count", "prior_p", "p_value"]]

    def _calculate_prior_probability(self):
        if self.config.pipeline_settings.use_hicpro_bias:
            self.data["prior_p"] = self.interChrProb * self.data["bias_1"] * self.data["bias_2"]
        else:
            self.data["prior_p"] = self.interChrProb

        if self.config.pipeline_settings.normalize_inter_expected_frequency:
            self.data["prior_p"] /= self.total_possible_interactions

    def _calculate_p_values(self):
        # Binomial distribution survival function for p-value calculation
        self.data["p_value"] = self.data.map_partitions(
            lambda df: binom.sf(
                df["interaction_count"] - 1,
                self.total_possible_interactions,
                df["prior_p"]
            ),
            meta=("p_value", "float64")
        )
