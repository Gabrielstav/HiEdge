# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from scipy.stats import nchypergeom_wallenius as nchg
from scipy.stats import binom
import dask.dataframe as dd
from scipy.interpolate import UnivariateSpline

# TODO:
#  later: implement nchg test for intra p-vals


class IntraPValueCalculator:
    def __init__(self, data: dd.DataFrame, spline: UnivariateSpline, metadata, config: Config):
        self.data = data
        self.spline = spline
        self.metadata = metadata
        self.config = config

    def run(self):
        print(f"Columns in data for p-val: {self.data.columns}")
        self._calculate_expected_frequency()
        self._calculate_p_values()
        return self.data

    def _calculate_expected_frequency(self):
        # Use the spline to calculate expected frequency based on genomic distance
        self.data["expected_frequency"] = self.data["genomic_distance"].map(lambda x: self.spline(x))

    def _calculate_p_values(self):
        # Adjust expected frequency with bias terms if configured
        if self.config.statistical_settings.use_hicpro_bias:
            self.data["adjusted_expected_frequency"] = self.data["expected_frequency"] * self.data["bias_1"] * self.data["bias_2"]
        else:
            self.data["adjusted_expected_frequency"] = self.data["expected_frequency"]

        total_interactions = self._recalculate_total_interactions()

        # Calculate p-values using a binomial test
        self.data["p_value"] = self.data.apply(lambda row: binom.sf(row["interaction_count"] - 1, total_interactions, row["adjusted_expected_frequency"]), axis=1)

    def _recalculate_total_interactions(self) -> int:
        # recalculate the total interactions present in the data (after filtering)
        print(self.metadata.interaction_count_per_chromosome_intra)
        total_possible_interactions = sum(self.metadata.interaction_count_per_chromosome_intra[chrom] for chrom in self.metadata.chromosomes_present)
        return total_possible_interactions

class InterPValueCalculator:

    def __init__(self, config: Config, data, metadata, total_possible_interactions):
        self.config = config
        self.data = data
        self.metadata = metadata
        self.interChrProb = data.interChrProb
        self.total_possible_interactions = total_possible_interactions

    def run(self):
        self._calculate_prior_probability()
        self._calculate_p_values()
        return self.data[["interaction_count", "prior_p", "p_value"]]

    def _calculate_prior_probability(self):
        if self.config.statistical_settings.use_hicpro_bias:
            self.data["prior_p"] = self.interChrProb * self.data["bias_1"] * self.data["bias_2"]
        else:
            self.data["prior_p"] = self.interChrProb

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
