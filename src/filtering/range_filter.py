# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
import dask.dataframe as dd
import pandas as pd

class FilterRanges:

    def __init__(self, config: Config):
        self.config = config
        self.region_type = self._determine_region_type()
        self.transformed_regions = self._transform_regions()
        self.cached_omit_dict = None
        self.cached_select_dict = None

    def _determine_region_type(self):
        if self.config.pipeline_settings.select_regions:
            return "select"
        elif self.config.pipeline_settings.omit_regions:
            return "omit"
        else:
            return None

    def _transform_regions(self):
        regions = self.config.pipeline_settings.select_regions if self.region_type == "select" else self.config.pipeline_settings.omit_regions
        rows = [{"chromosome": chromosome, "start": genomic_range.start, "end": genomic_range.end}
                for chromosome, ranges in regions.items() for genomic_range in ranges]
        return pd.DataFrame(rows)

    @staticmethod
    def _preprocess_regions(transformed_regions):
        region_dict = {}
        for row in transformed_regions.itertuples(index=False):
            region_dict.setdefault(row.chromosome, []).append((row.start, row.end))
        return region_dict

    def _get_omit_dict(self):
        if self.cached_omit_dict is None:
            self.cached_omit_dict = self._preprocess_regions(self.transformed_regions)
        return self.cached_omit_dict

    def _get_select_dict(self):
        if self.cached_select_dict is None:
            self.cached_select_dict = self._preprocess_regions(self.transformed_regions)
        return self.cached_select_dict

    @staticmethod
    def _check_overlap_batch(df, region_dict):
        omit_series = pd.Series(False, index=df.index)
        for chromosome, regions in region_dict.items():
            for start, end in regions:
                overlap_cond = (
                        ((df["chr_1"] == chromosome) & (df["start_1"] <= end) & (df["end_1"] >= start)) |
                        ((df["chr_2"] == chromosome) & (df["start_2"] <= end) & (df["end_2"] >= start))
                )
                omit_series = omit_series | overlap_cond
        return ~omit_series

    def filter_omit_regions(self, data: dd.DataFrame) -> dd.DataFrame:
        omit_dict = self._get_omit_dict()
        mask = data.map_partitions(
            lambda df: self._check_overlap_batch(df, omit_dict),
            meta=bool
        )
        return data[mask]

    def filter_select_regions(self, data: dd.DataFrame) -> dd.DataFrame:
        select_dict = self._get_select_dict()
        mask = data.map_partitions(
            lambda df: self._check_overlap_batch(df, select_dict),
            meta=bool
        )
        return data[mask]
