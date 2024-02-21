# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.config_loader import GenomicRange as Gr
import dask.dataframe as dd
import pandas as pd
from typing import Dict, List

class FilterRanges:

    def __init__(self, config: Config):
        self.config = config
        self.select_regions = config.pipeline_settings.select_regions
        self.omit_regions = config.pipeline_settings.omit_regions

    @staticmethod
    def _check_overlap(df: pd.DataFrame, regions: Dict[str, List[Gr]]):
        """Check if rows in df strictly overlap with any region in regions dict."""
        mask = pd.Series(False, index=df.index)
        for chrom, grs in regions.items():
            for gr in grs:
                # Adjusted overlap condition to exclude interactions where only the end matches the start of a range
                overlap_cond = (
                    ((df["chr_1"] == chrom) & (df["start_1"] < gr.end) & (df["end_1"] > gr.start)) |
                    ((df["chr_2"] == chrom) & (df["start_2"] < gr.end) & (df["end_2"] > gr.start))
                )
                mask |= overlap_cond
        return mask

    def filter_regions(self, data: dd.DataFrame, regions: Dict[str, List[Gr]], keep: bool = True) -> dd.DataFrame:
        """Filters data based on regions. If keep is True, keeps rows overlapping with regions; otherwise removes them."""

        if not isinstance(data, dd.DataFrame):
            data = dd.from_pandas(data, npartitions=1)

        mask = data.map_partitions(
            lambda partition: self._check_overlap(partition, regions),
            meta=bool
        )
        if keep:
            return data[mask]
        else:
            return data[~mask]

    def apply_filters(self, data: dd.DataFrame) -> dd.DataFrame:
        if self.config.pipeline_settings.select_specific_regions:
            # Keep rows that overlap with select regions
            return self.filter_regions(data, self.select_regions, keep=True)
        elif self.config.pipeline_settings.omit_specific_regions:
            # Remove rows that overlap with omit regions
            return self.filter_regions(data, self.omit_regions, keep=False)
        return data
