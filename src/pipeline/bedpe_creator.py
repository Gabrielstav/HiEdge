# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import BedpeOutput, GroupedFiles
import dask.dataframe as dd

# TODO: Validation of input data and metadata

class BedpeCreator:

    def __init__(self, config: Config, grouped_files: GroupedFiles):
        self.config = config
        self.grouped_files = grouped_files

    def _make_bedpe_dask(self) -> BedpeOutput:

        bed_dtypes = {'chr': str, 'start': 'int32', 'end': 'int32', 'idx': 'int32'}
        matrix_dtypes = {'id1': 'int32', 'id2': 'int32', 'interaction_count': 'int32'}

        bed_ddf = dd.read_csv(self.grouped_files.bed_file, sep="\t", header=None, names=["chr", "start", "end", "idx"], dtype=bed_dtypes)
        matrix_ddf = dd.read_csv(self.grouped_files.matrix_file, sep="\t", header=None, names=["id1", "id2", "interaction_count"], dtype=matrix_dtypes)

        merged_ddf = matrix_ddf.merge(bed_ddf, left_on="id1", right_on="idx", how="left") \
            .merge(bed_ddf, left_on="id2", right_on="idx", how="left", suffixes=("_1", "_2"))

        bedpe_ddf = merged_ddf[["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count"]]

        # Now create the BedpeOutput dataclass instance using the metadata and the bedpe_ddf
        bedpe_output = BedpeOutput(metadata=self.grouped_files.metadata, data=bedpe_ddf)
        return bedpe_output


