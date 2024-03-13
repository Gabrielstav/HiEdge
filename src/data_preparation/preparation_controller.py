# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.data_preparation.bias_merger import BiasMerger
from src.data_preparation.interaction_creator import InteractionCreator
from src.setup.containers import InteractionOutput
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

        return InteractionOutput(metadata=self.grouped_files.metadata, data=interaction_df)
