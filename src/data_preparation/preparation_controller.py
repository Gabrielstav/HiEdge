# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.data_preparation.interaction_calculator import InteractionCalculator
from src.data_preparation.bias_merger import BiasMerger
from src.data_preparation.interaction_creator import InteractionCreator
from src.setup.data_structures import InteractionOutput
from src.setup.config_loader import Config

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

    def calculate_interactions_based_on_config(self, interaction_df):
        if self.config.statistical_settings.normalize_intra_expected_frequency or self.config.statistical_settings.normalize_inter_expected_frequency:
            interaction_calculator = InteractionCalculator(self.config, interaction_df)
            self._calculate_and_update_interactions(interaction_calculator)

    def _calculate_and_update_interactions(self, interaction_calculator):
        # Initialize placeholders for results
        results = {
            "total_possible_intra": None,
            "total_possible_inter": None,
            "interaction_count_per_chromosome_intra": None,
            "interaction_count_per_chromosome_inter": None
        }

        if self.config.statistical_settings.normalize_intra_expected_frequency:
            results["total_possible_intra"], results["interaction_count_per_chromosome_intra"] = interaction_calculator.calculate_intra_interactions()

        if self.config.statistical_settings.normalize_inter_expected_frequency:
            results["total_possible_inter"], results["interaction_count_per_chromosome_inter"] = interaction_calculator.calculate_inter_interactions()

        self._update_metadata_with_interactions(results)

    def _update_metadata_with_interactions(self, results):
        # Directly update the metadata within grouped_files
        self.grouped_files.metadata.max_possible_interaction_count_intra = results["total_possible_intra"]
        self.grouped_files.metadata.max_possible_interaction_count_inter = results["total_possible_inter"]
        self.grouped_files.metadata.interaction_count_per_chromosome_intra = results.get("interaction_count_per_chromosome_intra")
        self.grouped_files.metadata.interaction_count_per_chromosome_inter = results.get("interaction_count_per_chromosome_inter")
