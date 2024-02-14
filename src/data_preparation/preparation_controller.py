# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.data_preparation.interaction_calculator import InteractionCalculator
from src.data_preparation.bias_merger import BiasMerger
from src.data_preparation.interaction_creator import InteractionCreator
from src.setup.data_structures import InteractionOutput

class DataPreparationController:

    def __init__(self, config, grouped_files):
        self.config = config
        self.grouped_files = grouped_files
        self.interaction_processor = InteractionCreator(config, grouped_files, bed_length=None)
        self.interaction_calculator = InteractionCalculator(config, grouped_files.matrix_file)

    def run(self):
        result = self.interaction_processor.create_interaction_dataframe()
        interaction_df = result[0]
        bed_length = result[1] if len(result) > 1 else None  # Unpack bed_length if available

        if self.config.statistical_settings.use_hicpro_bias:
            bias_merger = BiasMerger(self.config, interaction_df, self.grouped_files.metadata.bias_file_path, bed_length)
            interaction_df = bias_merger.process_and_merge_bias()

        self.interaction_calculator = InteractionCalculator(self.config, interaction_df)

        total_interactions = self.interaction_calculator.calculate_total_interactions()

        updated_metadata = self._update_metadata_with_interactions(total_interactions)
        return InteractionOutput(metadata=updated_metadata, data=interaction_df)

    def _update_metadata_with_interactions(self, total_interactions):
        # unpack interactions
        total_intra, total_inter, interaction_count_per_chromosome = total_interactions

        # update the metadata field in grouped_files
        self.grouped_files.metadata.total_interaction_count_intra = total_intra
        self.grouped_files.metadata.total_interaction_count_inter = total_inter
        self.grouped_files.metadata.interaction_count_per_chromosome = interaction_count_per_chromosome

        return self.grouped_files.metadata
