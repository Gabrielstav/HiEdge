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

        self.data = self.data.set_index("chr_1", sorted=True)
        # self._interaction_count_per_chromosome()
        self._total_count()
        self._calculate_expected_frequency()
        self._calculate_p_values()
        self._plot_p_values_vs_distance()
        print(f"Tail of data after p-values: {self.data.head(30)}")
        return self.data

    def _calculate_p_values(self):
        # Instead of direct assignment, use apply to calculate p-values row-wise
        self.data["p_value"] = self.data.apply(self._calculate_p_value, axis=1, meta=("p_value", float))
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

    def _plot_p_values_vs_distance(self):
        # Ensure the DataFrame is computed; consider using .compute() carefully if the dataset is large
        # You might want to sample the data instead if it's very large
        sample_data = self.data.sample(frac=1).compute() if isinstance(self.data, dd.DataFrame) else self.data

        plt.figure(figsize=(10, 6))
        plt.scatter(sample_data["genomic_distance"], sample_data["p_value"], alpha=0.5)
        plt.title('P-Values vs Genomic Distance')
        plt.xlabel('Genomic Distance')
        plt.ylabel('P-Value')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True, which="both", ls="--")
        plt.show()

class PlaceHolder:
    def _unique_bins_per_chromosome(self) -> dd.DataFrame:
        bins_per_chrom_series = pd.Series(self.unique_bins_per_chromosome, name="unique_bins_per_chromosome")
        bins_per_chrom_dd = dd.from_pandas(bins_per_chrom_series.reset_index().rename(columns={"index": "chr_1"}), npartitions=1)
        data = dd.merge(self.data, bins_per_chrom_dd, on="chr_1", how="left")
        return data

    def _calculate_p_values(self):
        # Vectorized calculation of p-values using binomial survival function
        # Ensure necessary columns like 'total_interactions' are available

        adjusted_expected_frequency = self.data["expected_frequency"]
        if self.config.statistical_settings.use_hicpro_bias:
            adjusted_expected_frequency *= self.data["bias_1"] * self.data["bias_2"]

    def something(self):
        adjusted_expected_frequency = self.data["expected_frequency"]
        if self.config.statistical_settings.use_hicpro_bias:
            adjusted_expected_frequency *= self.data["bias_1"] * self.data["bias_2"]



        self.data["p_value"] = binom.sf(
            self.data["interaction_count"] - 1,
            self.data["interactions_per_chromosome"],
            adjusted_expected_frequency
        )


class OldIntraPValueCalculator:
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
