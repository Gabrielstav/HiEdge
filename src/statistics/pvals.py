# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from scipy.stats import nchypergeom_wallenius as nchg
from scipy.stats import binom
import dask.dataframe as dd
import numpy as np
from dask.dataframe.core import Series
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import pandas as pd
pd.options.mode.chained_assignment = None

# TODO:
#  later: implement nchg test for intra p-vals

class IntraPValueCalculator:

    def __init__(self, data: dd.DataFrame, spline, metadata, config: Config, unique_bins_per_chromosome: dict):
        self.data = data
        self.spline = spline
        self.metadata = metadata
        self.config = config
        self.unique_bins_per_chromosome = unique_bins_per_chromosome

    def run(self):
        if not isinstance(self.data, dd.DataFrame):
            print("Data is not a Dask DataFrame, converting.")
            self.data = dd.from_pandas(self.data, npartitions=1)

        self.data = self.data.set_index("chr_1", sorted=True, drop=False)
        # self._interaction_count_per_chromosome()
        self._total_count()
        self._calculate_expected_frequency()
        self._calculate_p_values()
        # self._plot_p_values_vs_distance()
        print(f"Tail of data after p-values: {self.data.head(30)}")
        return self.data

    def _calculate_p_values(self):
        # Instead of direct assignment, use apply to calculate p-values row-wise
        self.data["p_value"] = self.data.apply(self._calculate_p_value, axis=1, meta=("p_value", float))
        # TODO: Check expected frequency, bias, interaction count, total interactions and adjusted expected frequency

        # number of successes, number of trials, and success probability in each trial
        # so for each row, the number of successes is the interaction count, the number of trials is the total interactions, and the success probability is the expected frequency
        # it's the y prob from the spline for a given distance, but something is wrong with this calculation since p-vals are too small...

    def _calculate_p_value(self, row):
        # Ensure inputs are scalars or numpy arrays
        try:
            # Calculate adjusted expected frequency within the function
            adjusted_expected_frequency = row["expected_frequency"]
            if self.config.statistical_settings.use_hicpro_bias:
                adjusted_expected_frequency *= row["bias_1"] * row["bias_2"]

            return binom.sf(row["interaction_count"] - 1, row["total_interactions"], adjusted_expected_frequency)

        except Exception as e:
            # Handle or log the exception
            print(f"Error calculating p-value: {e}")
            return np.nan  # Return a NaN or some sentinel value on error

    def _calculate_expected_frequency(self):
        # Utilize map_partitions to apply the spline's predict function
        # This ensures the operation remains lazy and is applied to each partition
        self.data["expected_frequency"] = self.data.map_partitions(
            lambda df: df.apply(lambda row: self.spline.predict(row["genomic_distance"]), axis=1),
            meta=("expected_frequency", np.float64)
        )

    def _interaction_count_per_chromosome(self) -> dd.DataFrame:
        contact_count = self.data.groupby("chr_1")["interaction_count"].transform("sum")
        self.data["interactions_per_chromosome"] = contact_count
        return self.data

    def _total_count(self):
        total_interactions = self.data["interaction_count"].sum()
        self.data["total_interactions"] = total_interactions
        return self.data


class InterPValueCalculator:
    def __init__(self, data: dd.DataFrame, config: Config, metadata, total_possible_interactions):
        self.config = config
        self.data = data
        self.metadata = metadata
        self.total_possible_interactions = total_possible_interactions

        # Verify interChrProb exists
        if "interChrProb" not in self.data.columns:
            raise ValueError("Input data missing interChrProb column")

    def run(self):
        self._calculate_prior_probability()
        self._calculate_p_values()
        return self.data

    def _calculate_prior_probability(self):
        if self.config.statistical_settings.use_hicpro_bias:
            self.data["prior_p"] = self.data["interChrProb"] * self.data["bias_1"] * self.data["bias_2"]
        else:
            self.data["prior_p"] = self.data["interChrProb"]

    def _calculate_p_values(self):
        # Use total observed sum as number of trials
        total_observed = self.data["totalObservedSum"].loc[0]

        self.data["p_value"] = self.data.map_partitions(
            lambda df: binom.sf(
                df["interaction_count"] - 1,
                total_observed,  # Changed from total_possible_interactions
                df["prior_p"]
            ),
            meta=("p_value", "float64")
        )
