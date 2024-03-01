# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
import dask.dataframe as dd
import numpy as np

class MidpointCalculator:

    def __init__(self, config: Config):
        self.config = config

    @staticmethod
    def calculate_midpoints(data: dd.DataFrame) -> dd.DataFrame:
        data["midpoint_1"] = ((data["start_1"] + data["end_1"]) / 2).astype("int32")
        data["midpoint_2"] = ((data["start_2"] + data["end_2"]) / 2).astype("int32")
        return data


class DistanceCalculator:

    def __init__(self, config: Config):
        self.config = config

    @staticmethod
    def calculate_distances(data: dd.DataFrame) -> dd.DataFrame:
        # Calculate the genomic distance
        data["genomic_distance"] = abs(data["midpoint_1"] - data["midpoint_2"]).astype("int32")

        # Filter out data with genomic distance == 0 (self-interactions)
        filtered_data = data[data["genomic_distance"] != 0]

        return filtered_data


class DistanceFilter:

    def __init__(self, config: Config, data: dd.DataFrame, resolution: int):
        self.config = config
        self.data = data
        self.resolution = resolution

    def filter_distances(self) -> dd.DataFrame:
        if self.resolution in self.config.pipeline_settings.interaction_distance_filters and self.config.pipeline_settings.use_interaction_distance_filters:
            distance_filter_config = self.config.pipeline_settings.interaction_distance_filters[self.resolution]
            lower_bound = distance_filter_config.min_distance
            upper_bound = distance_filter_config.max_distance
            self.data = self.data[(self.data["genomic_distance"] >= lower_bound) & (self.data["genomic_distance"] <= upper_bound)]
        return self.data


class EqualOccupancyBinner:

    def __init__(self, config: Config, input_data: dd.DataFrame):
        self.config = config
        self.input_data = input_data

    def bin_data(self, data: dd.DataFrame) -> dd.DataFrame:
        sorted_data = self._sort_data_by_distance(data)
        total_contacts = self._calculate_total_contacts(data)
        binned_data = self._assign_to_bins(sorted_data, total_contacts, self.config.statistical_settings.metabin_occupancy)
        bin_stats = self._calculate_metabin_stats(binned_data)

        print(f"Bin stats: {bin_stats}")
        print(f"Data columns: {data.columns}")
        return bin_stats

    @staticmethod
    def _sort_data_by_distance(input_data: dd.DataFrame) -> dd.DataFrame:
        return input_data.sort_values("genomic_distance")

    @staticmethod
    def _calculate_total_contacts(data: dd.DataFrame) -> int:
        return data["interaction_count"].sum()

    @staticmethod
    def _assign_to_bins(sorted_data: dd.DataFrame, total_contacts: int, num_bins: int) -> dd.DataFrame:
        target_per_bin = total_contacts / num_bins
        sorted_data["cumulative_count"] = sorted_data["interaction_count"].cumsum()
        sorted_data["metabin_label"] = (sorted_data["cumulative_count"] / target_per_bin).astype("int")

        # Handle the tiebreak condition
        sorted_data = EqualOccupancyBinner._apply_tiebreak_condition(sorted_data)

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
        # Directly update metabin_label with adjusted values
        data["metabin_label"] = np.where(condition, data["metabin_label"], data["metabin_label"].shift(1).fillna(-1).astype(int))

        data = data.drop(["needs_tiebreak"], axis=1)

        if not isinstance(data, dd.DataFrame):
            data = dd.from_pandas(data, npartitions=1)

        return data

    @staticmethod
    def _calculate_metabin_stats(data_binned: dd.DataFrame) -> dd.DataFrame:
        # Ensure changes are applied correctly by assigning the result back
        # data_binned = data_binned.drop(["bin_label"], axis=1, errors="ignore")

        bin_stats = data_binned.groupby("metabin_label").agg({
            "genomic_distance": "mean",
            "interaction_count": "sum"
        }).reset_index()

        bin_stats = bin_stats.rename(columns={
            "genomic_distance": "average_distance",
            "interaction_count": "total_contacts"
        })

        # calculate average contact probability per metabin
        bin_stats["average_contact_probability"] = bin_stats["total_contacts"] / bin_stats["total_contacts"].sum()

        return bin_stats

    @staticmethod
    def _merge_bin_stats_with_interaction_dataframe(bin_stats: dd.DataFrame, interaction_dataframe: dd.DataFrame) -> dd.DataFrame:

        # Merge bin statistics back into the original DataFrame (not used anymore)
        if not isinstance(bin_stats, dd.DataFrame):
            print("Converting bin_stats stats data to dask dataframe")
            bin_stats = dd.from_pandas(bin_stats, npartitions=1)
        if not isinstance(interaction_dataframe, dd.DataFrame):
            print("Converting interaction dataframe to dask dataframe")
            interaction_dataframe = dd.from_pandas(interaction_dataframe, npartitions=1)

        return dd.merge(interaction_dataframe, bin_stats, on="metabin_label", how="left")


class FrequencyAggregator:

    def __init__(self, config: Config):
        self.config = config

    # Aggregate frequencies of genomic distances
    @staticmethod
    def aggregate_frequencies(data: dd.DataFrame) -> dd.DataFrame:
        print(data)
        aggregated_data = data.groupby("average_distance").agg({
            "interaction_count": ["sum", "mean", "std"]
        }).reset_index()
        aggregated_data.columns = ["average_distance", "interaction_sum", "interaction_mean", "interaction_std"]
        return aggregated_data
