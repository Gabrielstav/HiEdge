# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import dask.dataframe as dd
from src.setup.config_loader import Config

class MidpointCalculator:

    def __init__(self, config: Config):
        self.config = config

    @staticmethod
    def calculate_midpoints(data: dd.DataFrame) -> dd.DataFrame:
        data["midpoint_1"] = ((data["start_1"] + data["end_1"]) / 2).astype("int32")
        data["midpoint_2"] = ((data["start_2"] + data["end_2"]) / 2).astype("int32")
        return data


class TotalInteractionsCalculator:

    def __init__(self, config: Config, data: dd.DataFrame):
        self.config = config
        self.data = data
        self.total_interactions: int = 0

    def calculate_total_interactions(self) -> int:
        total_interactions = self.data["interaction_count"].sum()
        self.total_interactions = total_interactions
        return total_interactions


class DistanceCalculator:

    def __init__(self, config: Config):
        self.config = config

    def calculate_distances(self, data: dd.DataFrame) -> dd.DataFrame:
        # Calculate the genomic distance
        data["genomic_distance"] = abs(data["midpoint_1"] - data["midpoint_2"]).astype("int32")

        # Filter out data with genomic distance == 0 (self-interactions)
        if self.config.pipeline_settings.filter_self_interactions:
            filtered_data = data[data["genomic_distance"] != 0]
        else:
            filtered_data = data

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
