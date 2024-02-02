# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
import dask.dataframe as dd


class DataPreparationController:

    def __init__(self, config: Config, input_data):
        self.config = config
        self.input_data = input_data

    def run(self):
        # Calculate inter stats
        inter_stats = InterStatisticsCalculator(self.input_data, self.config).compute_inter_stats()

        return inter_stats

class InterStatisticsCalculator:

    def __init__(self, data: dd.DataFrame, config: Config):
        self.data = data
        self.config = config

    def compute_inter_stats(self):
        # Filter for interchromosomal interactions
        inter_data = self.data[self.data["chr_1"] != self.data["chr_2"]]

        if self.data[self.data["chr_1"] == self.data["chr_2"]] in inter_data:
            raise ValueError("Inter dataset contains intrachromosomal interactions!!!")  # Just sanity check before adding validation

        # Count all interactions and compute interaction probability, add as column in data (or as metadata?)
        interdata_count = len(inter_data.index)
        self.data["interdata_count"] = self.calculate_interaction_probabilities(interdata_count)

        return interdata_count

    @staticmethod
    def calculate_interaction_probabilities(interdata_count):
        inter_probability = 1 / interdata_count
        return inter_probability

