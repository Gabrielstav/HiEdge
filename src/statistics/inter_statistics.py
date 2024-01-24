# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import FilteringOutput, BlacklistOutput, CytobandOutput
import dask.dataframe as dd

class ResolveBedpeInput:

    def __init__(self, config: Config, data_output):
        self.config = config
        self.data_output = data_output

    def resolve_input(self):
        # Mapping configuration options to data classes
        data_class_mapping = {
            "cytoband": (CytobandOutput, self.config.pipeline_settings.filter_cytobands),
            "blacklist": (BlacklistOutput, self.config.pipeline_settings.filter_blacklist)
        }

        for data_class, (output_type, use_filter) in data_class_mapping.items():
            if isinstance(self.data_output, output_type) and use_filter:
                return output_type

        # Return the data field of the resolved data class
        return self.data_output.data

class DataPreparationController:

    def __init__(self, config: Config, input_data):
        self.config = config
        self.input_data = input_data

    def run(self):
        # Resolve the input data class based on config
        resolved_data = ResolveBedpeInput(self.config, self.input_data).resolve_input()

        # Calculate inter stats
        inter_stats = InterStatisticsCalculator(resolved_data, self.config).compute_inter_stats()

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

