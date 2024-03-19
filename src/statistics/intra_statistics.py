# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.statistics.stat_utils import TotalInteractionsCalculator
import dask.dataframe as dd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class IntraStatsProcessor:

    def __init__(self, config, data: dd.DataFrame, unique_bins_per_chromosome: dict, resolution: int):
        self.config = config
        self.data = data
        self.total_interactions: int = self._calculate_total_interactions()
        self.unique_bins_per_chromosome = unique_bins_per_chromosome
        self.resolution = resolution

    def __post_init__(self):
        self.total_interactions = self._calculate_total_interactions()

    def run(self) -> dd.DataFrame:

        # metabin the data and calculate metabin statistics
        metabinned_data = MetaBinner(self.config, self.data, self.total_interactions).bin_data()

        binned_data_with_stats = self._calculate_metabin_stats(metabinned_data)

        # aggregate data to metabins
        binned_data_with_stats = binned_data_with_stats.groupby("metabin_label").agg({
            "unique_bins_in_metabin": "first",
            "total_possible_pairs_within_metabin": "first",
            "metabin_total_contact_count": "first",
            "average_contact_probability": "first",
            "average_genomic_distance": "first"
        }).reset_index()

        # sort data by average genomic distance
        binned_data_with_stats = binned_data_with_stats.sort_values("average_genomic_distance")

        # print binned sorted data
        print(f"binned sorted data: {binned_data_with_stats.head(5)}")

        return binned_data_with_stats

    def _calculate_metabin_stats(self, binned_data):
        metabin_stats_calculator = MetabinStatistics(input_data=binned_data, config=self.config, total_interaction_count=self.total_interactions, unique_bins_per_chromosome=self.unique_bins_per_chromosome, resolution=self.resolution)
        return metabin_stats_calculator.apply_metabin_statistics()

    def _calculate_total_interactions(self) -> int:
        total_interactions = TotalInteractionsCalculator(self.config, self.data).calculate_total_interactions()
        return total_interactions

class MetaBinner:

    def __init__(self, config: Config, input_data: dd.DataFrame, total_interactions: int):
        self.config = config
        self.data = input_data
        self.total_interactions = total_interactions

    def bin_data(self) -> dd.DataFrame:
        self._sort_data_by_distance(self.data)
        self._assign_to_metabins(self.data, self.total_interactions, self.config.statistical_settings.metabin_occupancy)
        print(f"Data after binning: {self.data}")
        return self.data

    @staticmethod
    def _assign_to_metabins(sorted_data: dd.DataFrame, total_interactions: int, number_of_equal_occupancy_bins: int) -> dd.DataFrame:

        number_of_equal_occupancy_bins = number_of_equal_occupancy_bins-1

        target_per_metabin = total_interactions // number_of_equal_occupancy_bins
        sorted_data["cumulative_count"] = sorted_data["interaction_count"].cumsum()
        sorted_data["metabin_label"] = (sorted_data["cumulative_count"] / target_per_metabin).astype("int")


        # Handle the tiebreak condition
        sorted_data = MetaBinner._apply_tiebreak_condition(sorted_data)
        print(f"Head of sorted data in assign_to_bins: {sorted_data.head(20)}")

        return sorted_data

    @staticmethod
    def _apply_tiebreak_condition(data: dd.DataFrame) -> dd.DataFrame:
        # Determine if tiebreak is needed
        data["needs_tiebreak"] = (
                (data["genomic_distance"] == data["genomic_distance"].shift(1)) &
                (data["metabin_label"] != data["metabin_label"].shift(1))
        )

        # Adjust bin labels where tiebreak is needed
        condition = ~data["needs_tiebreak"]
        data["metabin_label"] = np.where(condition, data["metabin_label"], data["metabin_label"].shift(1).fillna(-1).astype(int))

        data = data.drop(["needs_tiebreak"], axis=1)

        return data

    @staticmethod
    def _sort_data_by_distance(input_data: dd.DataFrame) -> dd.DataFrame:
        return input_data.sort_values("genomic_distance")


class MetabinStatistics:

    def __init__(self, input_data: dd.DataFrame, config: Config, total_interaction_count: int, unique_bins_per_chromosome: dict, resolution: int):
        self.original_data = input_data
        self.config = config
        self.total_interaction_count = total_interaction_count
        self.distance_scaling_factor: int = 1000000
        self.unique_bins_per_chromosome = unique_bins_per_chromosome
        self.resolution = resolution

    def apply_metabin_statistics(self):
        """
        Sequentially apply metabin statistics to the metabinned data.
        """
        if not isinstance(self.original_data, dd.DataFrame):
            print(f"Data is not a dask dataframe, converting. Data type: {type(self.original_data)}")
            data = dd.from_pandas(self.original_data, npartitions=1)
        else:
            data = self.original_data

        data = self._create_unique_bin_id(data)
        data = self._get_unique_bins_per_metabin(data)
        data = self._possible_pairs_per_metabin(data)
        data = self._possible_pairs_in_metabin_range(data)
        data = self._unique_bins_per_chromosome(data)
        data = self._sum_of_contact_counts(data)
        data = self._average_contact_probability(data)
        data = self._average_genomic_distance(data)


        return data

    @staticmethod
    def _create_unique_bin_id(data: dd.DataFrame) -> dd.DataFrame:
        data["bin_id"] = data["chr_1"] + data["start_1"].astype(str)
        return data

    @staticmethod
    def _get_unique_bins_per_metabin(data: dd.DataFrame) -> dd.DataFrame:
        unique_bins_counts = data.groupby("metabin_label")["bin_id"].nunique().compute().reset_index()
        data["unique_bins_in_metabin"] = data["metabin_label"].map(unique_bins_counts.set_index("metabin_label")["bin_id"])
        return data

    @staticmethod
    def _possible_pairs_per_metabin(data: dd.DataFrame) -> dd.DataFrame:
        # this is by definition the possible pairs found in the genomic range of each metabin
        data["total_possible_pairs_within_metabin"] = (data["unique_bins_in_metabin"] * (data["unique_bins_in_metabin"] + 1)) // 2
        return data

    def _possible_pairs_in_metabin_range(self, data: dd.DataFrame) -> dd.DataFrame:
        # data["distance_range_in_metabin"] = data.groupby("metabin_label")["genomic_distance"].transform("max") - data.groupby("metabin_label")["genomic_distance"].transform("min")
        data["max_midpoint_2_metabin"] = data.groupby("metabin_label")["midpoint_2"].transform("max").compute()
        data["min_midpoint_1_metabin"] = data.groupby("metabin_label")["midpoint_1"].transform("min").compute()
        data["distance_range_in_metabin"] = data["max_midpoint_2_metabin"] - data["min_midpoint_1_metabin"]
        data["max_possible_pairs_within_metabin"] = (data["distance_range_in_metabin"] / self.resolution) * (data["distance_range_in_metabin"] / self.resolution + 1) / 2
        print(f"Heads of data in possible_pairs_in_metabin_range: {data.head(5)}")
        print(f"Tail of data in possible_pairs_in_metabin_range: {data.tail(5)}")

        return data

    def _unique_bins_per_chromosome(self, data: dd.DataFrame) -> dd.DataFrame:
        bins_per_chrom_series = pd.Series(self.unique_bins_per_chromosome, name="unique_bins_per_chromosome")
        bins_per_chrom_dd = dd.from_pandas(bins_per_chrom_series.reset_index().rename(columns={"index": "chr_1"}), npartitions=1)
        data = dd.merge(data, bins_per_chrom_dd, on="chr_1", how="left")
        return data

    @staticmethod
    def _sum_of_contact_counts(data: dd.DataFrame) -> dd.DataFrame:
        contact_count = data.groupby("metabin_label")["interaction_count"].transform("sum")
        data["metabin_total_contact_count"] = contact_count
        return data

    def _average_contact_probability(self, data: dd.DataFrame) -> dd.DataFrame:
        data["average_contact_probability"] = ((data["metabin_total_contact_count"] / data["max_possible_pairs_within_metabin"]) / self.total_interaction_count)
        return data

    def _average_genomic_distance(self, data: dd.DataFrame) -> dd.DataFrame:
        data["sum_distances"] = data.groupby("metabin_label")["genomic_distance"].transform("sum") / (self.distance_scaling_factor * data["unique_bins_per_chromosome"])
        data["average_genomic_distance"] = self.distance_scaling_factor*(data["sum_distances"]/data["unique_bins_in_metabin"])
        return data
