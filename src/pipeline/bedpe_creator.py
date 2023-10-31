# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from pathlib import Path
from src.setup.config_loader import ConfigMapper as Config
import pandas as pd
import dask.dataframe as dd

class BedpeCreator:

    def __init__(self, config: Config, matrix_file: Path, bed_file: Path):
        self.config = config
        self.bed_file = bed_file
        self.matrix_file = matrix_file
        # might need config later for filtering


    def make_bedpe_pandas(self, matrix_file: Path, bed_file: Path) -> pd.DataFrame:
        """
        Use Pandas to make BEDPE files
        :param matrix_file: Matrix file from HiC-Pro grouped by HicProInputFilePreparer in pipeline_input
        :param bed_file: Bed file from HiC-Pro grouped by HicProInputFilePreparer in pipeline_input
        :return: BEDPE dataframe made from matrix and bed files
        """

        bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chr", "start", "end", "idx"])
        matrix_df = pd.read_csv(matrix_file, sep="\t", header=None, names=["id1", "id2", "interaction_count"])

        merged_df = matrix_df.merge(bed_df, left_on="id1", right_on="idx", how="left") \
            .merge(bed_df, left_on="id2", right_on="idx", how="left", suffixes=("", "_2"))

        bedpe_df = merged_df[["chr_1", "start_1", "end", "chr_2", "start_2", "end_2", "interaction_count"]]

        return bedpe_df

    def make_bedpe_dask(self, matrix_file: Path, bed_file: Path) -> dd.DataFrame:
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
