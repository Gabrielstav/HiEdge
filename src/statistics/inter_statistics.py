# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.data_preparation.interaction_calculator import PossiblePairsCalculator
import dask.dataframe as dd


class DataPreparationController:

    def __init__(self, config: Config, input_data):
        self.config = config
        self.input_data = input_data

    def run(self):
        inter_stats = InterStatisticsCalculator(self.input_data, self.config).compute_inter_stats()
        return inter_stats

class InterStatisticsCalculator:
    def __init__(self, data: dd.DataFrame, config: Config):
        self.data = data
        self.config = config

        # Ensure data is Dask DataFrame
        if not isinstance(self.data, dd.DataFrame):
            self.data = dd.from_pandas(self.data, npartitions=1)

    def compute_inter_stats(self) -> dd.DataFrame:
        # Filter for interchromosomal interactions
        inter_data = self.data[self.data["chr_1"] != self.data["chr_2"]]

        # Calculate total observed interactions
        total_observed_count = len(inter_data)  # number of interacting pairs
        total_observed_sum = inter_data["interaction_count"].sum()  # sum of all interactions

        # Calculate probability (using their method)
        inter_probability = 1.0 / total_observed_count if total_observed_count > 0 else 0

        # Add these as properties of the DataFrame
        inter_data = inter_data.assign(
            interChrProb=inter_probability,
            totalObservedSum=total_observed_sum
        )

        return inter_data
