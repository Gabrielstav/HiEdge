# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from scipy.stats import nchypergeom_wallenius as nchg
from scipy.stats import binom
import dask.dataframe as dd
import numpy as np
from dask.dataframe.core import Series
from scipy.interpolate import UnivariateSpline
import pandas as pd
pd.options.mode.chained_assignment = None

# TODO:
#  later: implement nchg test for intra p-vals


class IntraPValueCalculator:
    def __init__(self, data: dd.DataFrame, spline, metadata, config: Config):
        self.data = data
        self.spline = spline
        self.metadata = metadata
        self.config = config
        self.total_interactions_per_chromosome = self._calculate_total_interactions_per_chromosome()

    def run(self):
        if not isinstance(self.data, dd.DataFrame):
            print("Converting to dask dataframe")
            self.data = dd.from_pandas(self.data, npartitions=1)
        if isinstance(self.data, dd.DataFrame):
            print(f"Data already in dask dataframe: {type(self.data)}")

        print(f"Data before exp freq: {self.data.head(5)}")
        self._calculate_expected_frequency()
        print(f"Data after exp freq: {self.data.head(5)}")
        self._calculate_p_values()
        # print head of data
        return self.data

    def _calculate_p_values(self):
        # Ensure data is prepared with total interactions
        self._prepare_data_with_total_interactions(self.data)

        # Adjust expected frequency with bias terms if configured
        adjusted_expected_frequency = self.data["expected_frequency"]
        if self.config.statistical_settings.use_hicpro_bias:
            adjusted_expected_frequency *= self.data["bias_1"] * self.data["bias_2"]

        # Vectorized p-value calculation
        self.data["p_value"] = binom.sf(
            self.data["interaction_count"] - 1,
            self.data["total_interactions"],
            adjusted_expected_frequency
        )

    def _prepare_data_with_total_interactions(self, data: dd.DataFrame):
        # TODO: Use unique_bins_per_chromosome to calculate total interactions per chromosome to vectorize this and keep it lazily evaluated insead of merging,
        #   this is observed interaction count (per chromosome) so I can just group by chromosome and sum interaction_count, so just keep in dask

        # like this:
        # data["total_interactions_per_chromosome"] = data.groupby["chr_1"](data["interaction_count"]).sum()

        # Convert total interactions per chromosome to a DataFrame
        total_interactions_df = pd.DataFrame.from_dict(self.total_interactions_per_chromosome, orient="index", columns=["total_interactions"])
        total_interactions_df["chr_1"] = total_interactions_df.index

        # Merge this DataFrame with the main data
        self.data = dd.merge(self.data, total_interactions_df, on="chr_1", how="left")

    def _calculate_total_interactions_per_chromosome(self):
        total_interactions_per_chromosome = self.data.groupby("chr_1")["interaction_count"].sum()
        return total_interactions_per_chromosome.to_dict()

    def _calculate_expected_frequency(self):
        # self.spline should be an instance of SplineFitter
        if isinstance(self.data, dd.DataFrame):
            self.data["expected_frequency"] = self.data["genomic_distance"].map(
                lambda x: self.spline.predict(x),
                meta=pd.Series(dtype=np.float64)
            )
        else:
            self.data["expected_frequency"] = self.data["genomic_distance"].apply(
                lambda x: self.spline.predict(x)
            )




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
