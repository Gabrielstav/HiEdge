# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import ConfigMapper as Config

# TODO: Pseudocode implementation, we need to figure out inter vs intra thing, and handle these separately

class PrepareInteractionData:
    def __init__(self, config: Config, data_input):
        self.config = config
        self.data_ddf = data_input.bedpe_ddf
        self.metadata = data_input.metadata

    def split_interactions(self):
        # Separate interactions into intra- and interchromosomal
        intra_df = self.data_ddf[self.data_ddf['chr_1'] == self.data_ddf['chr_2']]
        inter_df = self.data_ddf[self.data_ddf['chr_1'] != self.data_ddf['chr_2']]

        return intra_df, inter_df

    def calculate_genomic_distance(self, intra_df):
        # Calculate genomic distance only for intrachromosomal interactions
        intra_df['genomic_distance'] = abs(intra_df['midpoint_1'] - intra_df['midpoint_2'])
        return intra_df

    def process_data(self):
        intra_df, inter_df = self.split_interactions()
        intra_df = self.calculate_genomic_distance(intra_df)

        # Prepare StatInput dataclass for both types of interactions
        intra_stat_input = StatInput(metadata=self.metadata, bedpe_ddf=intra_df)
        inter_stat_input = StatInput(metadata=self.metadata, bedpe_ddf=inter_df)

        return intra_stat_input, inter_stat_input

