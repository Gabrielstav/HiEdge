# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import GroupedFiles
from src.data_preparation.validation import validate_interaction_data, validate_matrix_data, validate_bed_data
import dask.dataframe as dd

class InteractionCreator:

    def __init__(self, config: Config, grouped_files: GroupedFiles):
        self.config = config
        self.grouped_files = grouped_files

    def create_interaction_dataframe(self) -> dd.DataFrame:
        # Read the BED file into a Dask DataFrame with specified dtypes
        bed_ddf = dd.read_csv(
            self.grouped_files.bed_file,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "idx"]
        )

        # Read the matrix file into a Dask DataFrame with specified dtypes
        interaction_count_dtype = float if self.config.pipeline_settings.iced_data else int
        matrix_ddf = dd.read_csv(
            self.grouped_files.matrix_file,
            sep="\t",
            header=None,
            names=["id1", "id2", "interaction_count"],
            dtype={"id1": int, "id2": int, "interaction_count": interaction_count_dtype}
        )

        # Perform validation on bed and matrix dataframes
        validate_bed_data(bed_ddf)
        validate_matrix_data(matrix_ddf, self.config)

        # Merge matrix and bed dataframes to create the interaction dataframe
        merged_ddf = self._merge_dataframes(bed_ddf, matrix_ddf)

        # Perform validation on the interaction dataframe
        validate_interaction_data(merged_ddf, self.config)

        # Optionally round interaction count for iced data
        if self.config.pipeline_settings.iced_data and self.config.pipeline_settings.round_iced_matrices:
            merged_ddf["interaction_count"] = merged_ddf["interaction_count"].round().astype(int)

        return merged_ddf

    @staticmethod
    def _merge_dataframes(bed_ddf: dd.DataFrame, matrix_ddf: dd.DataFrame) -> dd.DataFrame:
        """Helper function to merge bed and matrix dataframes."""
        merged_ddf = matrix_ddf.merge(bed_ddf, left_on="id1", right_on="idx", how="left").rename(columns={"chr": "chr_1", "start": "start_1", "end": "end_1", "idx": "idx_1"})
        merged_ddf = merged_ddf.merge(bed_ddf, left_on="id2", right_on="idx", how="left").rename(columns={"chr": "chr_2", "start": "start_2", "end": "end_2", "idx": "idx_2"})
        return merged_ddf[["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count", "idx_1", "idx_2"]]



