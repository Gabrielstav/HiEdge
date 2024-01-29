# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import InteractionOutput, GroupedFiles
from src.data_preparation.interaction_calculator import InteractionCalculator
from src.data_preparation.bias_merger import BiasMerger
import dask.dataframe as dd


class InteractionCreator:

    def __init__(self, config: Config, grouped_files: GroupedFiles):
        self.config = config
        self.grouped_files = grouped_files

    def run(self):
        interaction_ddf = self.create_interaction_dataframe()

        if self.config.pipeline_settings.use_hicpro_bias:
            bias_series = BiasMerger.merge_bias(self.grouped_files, interaction_ddf)
            interaction_ddf = BiasMerger.merge_bias(interaction_ddf, bias_series)

        interaction_calculator = InteractionCalculator(self.config, interaction_ddf)
        total_intra, total_inter, interaction_count_per_chromosome = interaction_calculator.calculate_total_interactions()

        # evaluate delayed values
        total_intra, total_inter, interaction_count_per_chromosome = dd.compute(total_intra, total_inter, interaction_count_per_chromosome)

        # update metadata with interaction counts
        self.grouped_files.metadata.max_possible_interaction_count_intra = total_intra
        self.grouped_files.metadata.max_possible__interaction_count_inter = total_inter
        self.grouped_files.metadata.interaction_count_per_chromosome = interaction_count_per_chromosome

        return InteractionOutput(metadata=self.grouped_files.metadata, data=interaction_ddf)

    def create_interaction_dataframe(self) -> dd.DataFrame:

        bed_dtypes = {'chr': str, 'start': 'int32', 'end': 'int32', 'idx': 'int32'}
        matrix_dtypes = {'id1': 'int32', 'id2': 'int32', 'interaction_count': 'int32'}

        bed_ddf = dd.read_csv(self.grouped_files.bed_file, sep="\t", header=None, names=["chr", "start", "end", "idx"], dtype=bed_dtypes)
        matrix_ddf = dd.read_csv(self.grouped_files.matrix_file, sep="\t", header=None, names=["id1", "id2", "interaction_count"], dtype=matrix_dtypes)

        merged_ddf = matrix_ddf.merge(bed_ddf, left_on="id1", right_on="idx", how="left") \
            .merge(bed_ddf, left_on="id2", right_on="idx", how="left", suffixes=("_1", "_2"))

        bedpe_ddf = merged_ddf[["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count"]]

        return bedpe_ddf
