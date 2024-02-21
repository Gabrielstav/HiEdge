# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.data_preparation.interaction_calculator import InteractionCalculator
from src.data_preparation.bias_merger import BiasMerger
from src.data_preparation.interaction_creator import InteractionCreator
from src.setup.data_structures import InteractionOutput
from src.setup.config_loader import Config
from typing import Dict
import dask.dataframe as dd

class DataPreparationController:

    def __init__(self, config: Config, grouped_files):
        self.config = config
        self.grouped_files = grouped_files

    def run(self):
        # Create the interaction dataframe
        interaction_creator = InteractionCreator(self.config, self.grouped_files)
        interaction_df = interaction_creator.create_interaction_dataframe()

        # Merge bias with interaction data if necessary
        if self.config.statistical_settings.use_hicpro_bias:
            bias_merger = BiasMerger(self.config, interaction_df, self.grouped_files.metadata.bias_file_path, None)
            interaction_df = bias_merger.process_and_merge_bias()

        # Validate and calculate interactions based on config
        self.calculate_interactions_based_on_config(interaction_df)

        return InteractionOutput(metadata=self.grouped_files.metadata, data=interaction_df)

    def calculate_interactions_based_on_config(self, interaction_df: dd.DataFrame):
        if self.config.statistical_settings.normalize_expected_frequency:
            self._calculate_and_update_interactions(interaction_df)

    def _calculate_and_update_interactions(self, interaction_df: dd.DataFrame):
        # Pass the interaction dataframe to the InteractionCalculator at initialization
        interaction_calculator = InteractionCalculator(interaction_df)

        # Initialize placeholders for results
        total_intra = total_inter = interaction_count_per_chromosome_intra = None

        if self.config.pipeline_settings.interaction_type in ["mixed", "intra"]:
            # Now calling the method without passing the dataframe
            total_intra, interaction_count_per_chromosome_intra = interaction_calculator.calculate_intra_interactions()

        if self.config.pipeline_settings.interaction_type in ["mixed", "inter"]:
            # Now calling the method without passing the dataframe
            total_inter = interaction_calculator.calculate_inter_interactions()

        # Update metadata based on the calculated values
        self._update_metadata_with_interactions(total_intra, total_inter, interaction_count_per_chromosome_intra)

    def _update_metadata_with_interactions(self, total_intra: int, total_inter: int, interaction_count_per_chromosome_intra: Dict[str, int]):
        # Directly update the metadata within grouped_files
        self.grouped_files.metadata.max_possible_interaction_count_intra = total_intra
        self.grouped_files.metadata.max_possible_interaction_count_inter = total_inter
        self.grouped_files.metadata.interaction_count_per_chromosome_intra = interaction_count_per_chromosome_intra
