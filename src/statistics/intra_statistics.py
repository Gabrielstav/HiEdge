# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import SplineInput
import dask.dataframe as dd

class DataPreparationController:

    def __init__(self, config: Config, input_data):
        self.config = config
        self.input_data = input_data

    def run(self):
        # Prepare data for spline fitting
        prepared_data = self._prepare_data_for_spline(self.input_data)

        # Filter based on genomic distances
        filtered_data = self._filter_genomic_distances(prepared_data)

        # Trigger computation and instantiate SplineInput dataclass
        computed_data = filtered_data.compute() if isinstance(filtered_data, dd.DataFrame) else filtered_data
        return SplineInput(metadata=self.input_data.metadata, data=computed_data)

    def _filter_genomic_distances(self, data: dd.DataFrame) -> dd.DataFrame:
        resolution = self.input_data.metadata.resolution
        if resolution in self.config.pipeline_settings.interaction_distance_filters:
            distance_filter = DistanceFilter(self.config, data, resolution)
            return distance_filter.filter_distances()

        return data

    def _prepare_data_for_spline(self, data: dd.DataFrame) -> dd.DataFrame:
        midpoint_calculator = MidpointCalculator(self.config)
        data_with_midpoints = midpoint_calculator.calculate_midpoints(data)

        distance_calculator = DistanceCalculator(self.config)
        data_with_distances = distance_calculator.calculate_distances(data_with_midpoints)

        binner = EqualOccupancyBinner(self.config, data_with_distances)
        binned_data = binner.bin_data(data_with_distances)

        aggregator = FrequencyAggregator(self.config)
        aggregated_data = aggregator.aggregate_frequencies(binned_data)

        if self.input_data.metadata.resolution in self.config.pipeline_settings.interaction_distance_filters:
            distance_filter = DistanceFilter(self.config, aggregated_data, self.input_data.metadata.resolution)
            filtered_data = distance_filter.filter_distances()
            return filtered_data

        return aggregated_data


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
        data["genomic_distance"] = abs(data["midpoint_1"] - data["midpoint_2"]).astype("int32")
        return data


class EqualOccupancyBinner:

    def __init__(self, config, input_data: dd.DataFrame):
        self.config = config
        self.input_data = input_data

    @staticmethod
    def sort_data_by_distance(input_data: dd.DataFrame) -> dd.DataFrame:
        return input_data.sort_values("genomic_distance")

    @staticmethod
    def calculate_total_contacts(data: dd.DataFrame) -> int:
        return data["interaction_count"].sum()

    @staticmethod
    def assign_to_bins(sorted_data: dd.DataFrame, total_contacts: dd.DataFrame, num_bins: dd.DataFrame) -> dd.DataFrame:
        target_per_bin = total_contacts / num_bins
        sorted_data["cumulative_count"] = sorted_data["interaction_count"].cumsum()
        sorted_data["bin_label"] = (sorted_data["cumulative_count"] / target_per_bin).astype("int")

        # Handle the tiebreak condition
        sorted_data = EqualOccupancyBinner.apply_tiebreak_condition(sorted_data)

        return sorted_data

    @staticmethod
    def apply_tiebreak_condition(data: dd.DataFrame) -> dd.DataFrame:
        # Logic for applying the tiebreak condition
        data["needs_tiebreak"] = (
            (data["genomic_distance"] == data["genomic_distance"].shift(1)) &
            (data["bin_label"] != data["bin_label"].shift(1))
        )

        # Adjust bin labels where tiebreak is needed
        data["adjusted_bin_label"] = data["bin_label"].where(~data["needs_tiebreak"], data["bin_label"].shift(1))

        # Clean up and finalize
        data = data.drop(["needs_tiebreak"], axis=1)
        data = data.rename(columns={"adjusted_bin_label": "bin_label"})

        return data

    @staticmethod
    def calculate_bin_statistics(data_binned: dd.DataFrame) -> dd.DataFrame:
        # Calculate average genomic distance and contact probability for each bin
        bin_stats = data_binned.groupby("bin_label").agg({
            "genomic_distance": "mean",
            "interaction_count": ["mean", "sum"]
        }).compute()
        bin_stats.columns = ["average_distance", "average_contact_probability", "total_contacts"]
        return bin_stats

    def bin_data(self, data: dd.DataFrame) -> dd.DataFrame:
        sorted_data = self.sort_data_by_distance(data)
        total_contacts = self.calculate_total_contacts(data)
        binned_data = self.assign_to_bins(sorted_data, total_contacts, self.config.number_of_bins)
        bin_stats = self.calculate_bin_statistics(binned_data)
        return bin_stats

class FrequencyAggregator:

    def __init__(self, config: Config):
        self.config = config

    # Aggregate frequencies of genomic distances
    @staticmethod
    def aggregate_frequencies(data: dd.DataFrame) -> dd.DataFrame:
        aggregated_data = data.groupby("genomic_distance").agg({
            "interaction_count": ["sum", "mean", "std"]
        }).reset_index()
        aggregated_data.columns = ["genomic_distance", "interaction_sum", "interaction_mean", "interaction_std"]
        return aggregated_data


class DistanceFilter:

    def __init__(self, config: Config, data: dd.DataFrame, resolution: int):
        self.config = config
        self.data = data
        self.resolution = resolution

    def filter_distances(self) -> dd.DataFrame:
        if self.resolution in self.config.pipeline_settings.interaction_distance_filters and self.config.pipeline_settings.use_interaction_distance_filters:
            lower_bound = self.config.pipeline_settings.interaction_distance_filters[self.resolution]["lower"]
            upper_bound = self.config.pipeline_settings.interaction_distance_filters[self.resolution]["upper"]
            self.data = self.data[(self.data["genomic_distance"] >= lower_bound) & (self.data["genomic_distance"] <= upper_bound)]

        return self.data
