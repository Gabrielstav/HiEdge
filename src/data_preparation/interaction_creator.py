# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import InteractionOutput, GroupedFiles
from src.data_preparation.interaction_calculator import InteractionCalculator
from src.data_preparation.bias_merger import BiasMerger
import dask.dataframe as dd
from typing import Any


class InteractionCreator:

    def __init__(self, config: Config, grouped_files: GroupedFiles, bed_length):
        self.config = config
        self.grouped_files = grouped_files
        self.bed_length = bed_length

    def run(self):
        interaction_ddf = self.create_interaction_dataframe()
        print(interaction_ddf.head())  # TODO: Remove after debugging

        if self.config.pipeline_settings.use_hicpro_bias:
            bias_series = BiasMerger.merge_bias_with_interaction(self.grouped_files, interaction_ddf, self.bed_length)
            interaction_ddf = BiasMerger.merge_bias_with_interaction(interaction_ddf, bias_series, self.bed_length)

        interaction_calculator = InteractionCalculator(self.config, interaction_ddf)
        total_intra, total_inter, interaction_count_per_chromosome = interaction_calculator.calculate_total_interactions()

        # evaluate delayed values
        total_intra, total_inter, interaction_count_per_chromosome = dd.compute(total_intra, total_inter, interaction_count_per_chromosome)

        # update metadata with interaction counts
        self.grouped_files.metadata.max_possible_interaction_count_intra = total_intra
        self.grouped_files.metadata.max_possible_interaction_count_inter = total_inter
        self.grouped_files.metadata.interaction_count_per_chromosome = interaction_count_per_chromosome

        return InteractionOutput(metadata=self.grouped_files.metadata, data=interaction_ddf)

    def create_interaction_dataframe(self) -> Any:
        bed_dtypes = {"chr": str, "start": int, "end": int, "idx": int}
        matrix_dtypes = {"id1": int, "id2": int, "interaction_count": int}

        bed_ddf = dd.read_csv(self.grouped_files.bed_file, sep="\t", header=None, names=["chr", "start", "end", "idx"], dtype=bed_dtypes, encoding="latin-1")
        matrix_ddf = dd.read_csv(self.grouped_files.matrix_file, sep="\t", header=None, names=["id1", "id2", "interaction_count"], dtype=matrix_dtypes, encoding="latin-1")

        # Merge and explicitly rename columns to include _1 and _2 suffixes
        merged_ddf = matrix_ddf.merge(bed_ddf, left_on="id1", right_on="idx", how="left").rename(columns={"chr": "chr_1", "start": "start_1", "end": "end_1", "idx": "idx_1"})
        merged_ddf = merged_ddf.merge(bed_ddf, left_on="id2", right_on="idx", how="left").rename(columns={"chr": "chr_2", "start": "start_2", "end": "end_2", "idx": "idx_2"})

        # Now selecting the columns should work as expected
        bedpe_ddf = merged_ddf[["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count", "idx_1", "idx_2"]]

        # Conditionally calculate bed_length
        if self.config.statistical_settings.use_hicpro_bias:
            bed_length = bed_ddf.index.size.compute()
            return bedpe_ddf, bed_length
        else:
            return bedpe_ddf, None  # Return None for bed_length if bias data is not used
