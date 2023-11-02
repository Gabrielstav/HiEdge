# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from pathlib import Path
from src.setup.config_loader import ConfigMapper as Config
import pandas as pd
import dask.dataframe as dd

# TODO:
#   Infer metadata for each BEDPE file, save as headers in dfs?
#   Need to save: Resolution, normalization status, interaction type etc
#   Or is this not needed as each input file matching is parsed to one BEDPE DF, so the data is grouped correctly either way.
#   We only need to set window size, but that can be inferred from head of DF. But we also need to name output files,
#   and having metadata tags seems nice as this allows for expandability and easier testing.

class BedpeCreator:

    def __init__(self, config: Config, matrix_file: Path, bed_file: Path):
        self.config = config
        self.bed_file = bed_file
        self.matrix_file = matrix_file
        # might need config later for filtering

    @staticmethod
    def make_bedpe_dask(matrix_file: Path, bed_file: Path) -> dd.DataFrame:
        """
        Use Dask to make BEDPE files
        :param matrix_file: Matrix file from HiC-Pro
        :param bed_file: Bed file from HiC-Pro
        :return: BEDPE dataframe made from matrix and bed files
        """

        bed_ddf = dd.read_csv(bed_file, sep="\t", header=None, names=["chr", "start", "end", "idx"])
        matrix_ddf = dd.read_csv(matrix_file, sep="\t", header=None, names=["id1", "id2", "interaction_count"])

        merged_ddf = matrix_ddf.merge(bed_ddf, left_on="id1", right_on="idx", how="left") \
            .merge(bed_ddf, left_on="id2", right_on="idx", how="left", suffixes=("_1", "_2"))

        bedpe_ddf = merged_ddf[["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count"]]
        return bedpe_ddf


