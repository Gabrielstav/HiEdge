# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.filtering.bedpe_creator import BedpeOutput
from src.setup.data_structures import Metadata, FilteringOutput
from dask import dataframe as dd
import pandas as pd

class SplitByInteractionType:

    def __init__(self, bedpe_ddf: dd.DataFrame):
        self.bedpe_ddf = bedpe_ddf

    def split_datasets_by_interaction_type(self):
        intra_df = self.bedpe_ddf[self.bedpe_ddf["chr_1"] == self.bedpe_ddf["chr_2"]]
        inter_df = self.bedpe_ddf[self.bedpe_ddf["chr_1"] != self.bedpe_ddf["chr_2"]]
        return intra_df, inter_df

class SplitByInteractionType_old:

    def __init__(self, config: Config, bedpe_output: BedpeOutput):
        self.config = config
        self.bedpe_metadata = bedpe_output.metadata
        self.bedpe_ddf = bedpe_output.data

    def split_datasets_by_interaction_type(self):
        intra_df = self.bedpe_ddf[self.bedpe_ddf["chr_1"] == self.bedpe_ddf["chr_2"]]
        inter_df = self.bedpe_ddf[self.bedpe_ddf["chr_1"] != self.bedpe_ddf["chr_2"]]

        intra_metadata = Metadata(experiment=self.bedpe_metadata.experiment,
                                  resolution=self.bedpe_metadata.resolution,
                                  interaction_type="intra",
                                  bias_file_path=self.bedpe_metadata.bias_file_path)

        inter_metadata = Metadata(experiment=self.bedpe_metadata.experiment,
                                  resolution=self.bedpe_metadata.resolution,
                                  interaction_type="inter",
                                  bias_file_path=self.bedpe_metadata.bias_file_path)

        intra_output = FilteringOutput(data=intra_df, metadata=intra_metadata)
        inter_output = FilteringOutput(data=inter_df, metadata=inter_metadata)

        return intra_output, inter_output


class FilterChromosomes:

    def __init__(self, config: Config):
        self.config = config
        self.chromosome_action = self._determine_chromosome_action()
        self.chromosomes = self._get_chromosomes()

    def _determine_chromosome_action(self):
        if self.config.pipeline_settings.select_chromosomes:
            return "select"
        elif self.config.pipeline_settings.remove_chromosomes:
            return "remove"
        else:
            return None

    def _get_chromosomes(self):
        return self.config.pipeline_settings.select_chromosomes if self.chromosome_action == "select" \
            else self.config.pipeline_settings.remove_chromosomes

    def filter_data(self, data: dd.DataFrame) -> dd.DataFrame:
        if self.chromosome_action == "select":
            return self._filter_select_chromosomes(data)
        elif self.chromosome_action == "remove":
            return self._filter_remove_chromosomes(data)
        return data

    def _filter_select_chromosomes(self, data: dd.DataFrame) -> dd.DataFrame:
        is_intra = data["interaction_type"] == "intra"
        is_inter = data["interaction_type"] == "inter"
        in_chromosomes = data["chr_1"].isin(self.chromosomes) | data["chr_2"].isin(self.chromosomes)

        return data[(is_intra & (data["chr_1"].isin(self.chromosomes))) |
                    (is_inter & in_chromosomes)]

    def _filter_remove_chromosomes(self, data: dd.DataFrame) -> dd.DataFrame:
        is_intra = data["interaction_type"] == "intra"
        is_inter = data["interaction_type"] == "inter"
        not_in_chromosomes = ~data["chr_1"].isin(self.chromosomes) & ~data["chr_2"].isin(self.chromosomes)

        return data[(is_intra & (~data["chr_1"].isin(self.chromosomes))) |
                    (is_inter & not_in_chromosomes)]


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

class FilterInteractionDistances:

    def __init__(self, config: Config):
        self.config = config
        self.min_distance = self.config.pipeline_settings.min_interaction_range
        self.max_distance = self.config.pipeline_settings.max_interaction_range

    def filter_data(self, data: dd.DataFrame) -> dd.DataFrame:
        return data[(data["genomic_distance"] >= self.min_distance) & (data["genomic_distance"] <= self.max_distance)]


class FilterBias:
    @staticmethod
    def filter_bias(data: dd.DataFrame) -> dd.DataFrame:
        # Filter out interactions where either bias value is -1
        return data[(data["bias_1"] != -1) & (data["bias_2"] != -1)]





