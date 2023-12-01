# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.pipeline.bedpe_creator import BedpeOutput
from src.setup.data_structures import Metadata, SplittingOutput, FilteringOutput
from dataclasses import dataclass
from dask import dataframe as dd
import pandas as pd


class FilteringController:

    def __init__(self, config: Config, data: dd.DataFrame, metadata: Metadata):
        self.config = config
        self.data = data
        self.metadata = metadata

    def run_filters(self):
        # Split data based on interaction type
        splitter = SplitInteractions(self.config, BedpeOutput(data=self.data, metadata=self.metadata))
        intra_output, inter_output = splitter.split_datasets_by_interaction_type()

        # Apply chromosome filtering to both intra and inter datasets if specified
        if self.config.pipeline_settings.remove_chromosomes or self.config.pipeline_settings.select_chromosomes:
            chromosome_filter = FilterChromosomes(self.config)
            intra_output.data = chromosome_filter.filter_data(intra_output.data)
            inter_output.data = chromosome_filter.filter_data(inter_output.data)

        # Apply range filtering to both intra and inter datasets if specified
        range_filter = FilterRanges(self.config)
        if self.config.pipeline_settings.omit_regions:
            intra_output.data = range_filter.filter_omit_regions(intra_output.data)
            inter_output.data = range_filter.filter_omit_regions(inter_output.data)
        elif self.config.pipeline_settings.select_regions:
            intra_output.data = range_filter.filter_select_regions(intra_output.data)
            inter_output.data = range_filter.filter_select_regions(inter_output.data)

        # Return the filtered intra and inter datasets in FilteringOutput dataclasses
        return FilteringOutput(data=intra_output.data, metadata=intra_output.metadata), \
            FilteringOutput(data=inter_output.data, metadata=inter_output.metadata)


class SplitInteractions:

    def __init__(self, config: Config, bedpe_output: BedpeOutput):
        self.config = config
        self.bedpe_metadata = bedpe_output.metadata
        self.bedpe_ddf = bedpe_output.data

    def split_datasets_by_interaction_type(self):
        """
        Splits the bedpe dask dataframes by interaction type (inter- and intrachromosomal) to separate dask dataframes
        Each instance of the bedpe output dataclass contains one dask dataframe that potentially encompasses both inter/intra
        :return: Two dask dataframes, one containing only interchromosomal interaction and the other only intrachromosomal interactions
        """

        intra_df = self.bedpe_ddf[self.bedpe_ddf["chr_1"] == self.bedpe_ddf["chr_2"]]
        inter_df = self.bedpe_ddf[self.bedpe_ddf["chr_1"] != self.bedpe_ddf["chr_2"]]

        # Create metadata for intrachromosomal interactions
        intra_metadata = Metadata(
            experiment=self.bedpe_metadata.experiment,
            resolution=self.bedpe_metadata.resolution,
            interaction_type='intra'
        )

        # Create metadata for interchromosomal interactions
        inter_metadata = Metadata(
            experiment=self.bedpe_metadata.experiment,
            resolution=self.bedpe_metadata.resolution,
            interaction_type='inter'
        )

        # Create FilteringOutput instances
        intra_output = SplittingOutput(data=intra_df, metadata=intra_metadata)
        inter_output = SplittingOutput(data=inter_df, metadata=inter_metadata)

        return intra_output, inter_output


class FilterChromosomes:

    def __init__(self, config: Config):
        self.config = config
        self.chromosome_action = self._determine_chromosome_action()
        self.chromosomes = self._get_chromosomes()

    def _determine_chromosome_action(self):
        if self.config.pipeline_settings.select_chromosomes:
            return 'select'
        elif self.config.pipeline_settings.remove_chromosomes:
            return 'remove'
        else:
            return None

    def _get_chromosomes(self):
        return self.config.pipeline_settings.select_chromosomes if self.chromosome_action == 'select' \
            else self.config.pipeline_settings.remove_chromosomes

    def filter_data(self, data: dd.DataFrame) -> dd.DataFrame:
        if self.chromosome_action == 'select':
            return self._filter_select_chromosomes(data)
        elif self.chromosome_action == 'remove':
            return self._filter_remove_chromosomes(data)
        return data

    def _filter_select_chromosomes(self, data: dd.DataFrame) -> dd.DataFrame:
        is_intra = data['interaction_type'] == 'intra'
        is_inter = data['interaction_type'] == 'inter'
        in_chromosomes = data['chr_1'].isin(self.chromosomes) | data['chr_2'].isin(self.chromosomes)

        return data[(is_intra & (data['chr_1'].isin(self.chromosomes))) |
                    (is_inter & in_chromosomes)]

    def _filter_remove_chromosomes(self, data: dd.DataFrame) -> dd.DataFrame:
        is_intra = data['interaction_type'] == 'intra'
        is_inter = data['interaction_type'] == 'inter'
        not_in_chromosomes = ~data['chr_1'].isin(self.chromosomes) & ~data['chr_2'].isin(self.chromosomes)

        return data[(is_intra & (~data['chr_1'].isin(self.chromosomes))) |
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
            return 'select'
        elif self.config.pipeline_settings.omit_regions:
            return 'omit'
        else:
            return None

    def _transform_regions(self):
        regions = self.config.pipeline_settings.select_regions if self.region_type == 'select' else self.config.pipeline_settings.omit_regions
        rows = [{'chromosome': chromosome, 'start': genomic_range.start, 'end': genomic_range.end}
                for chromosome, ranges in regions.items() for genomic_range in ranges]
        return pd.DataFrame(rows)

    def _preprocess_regions(self, transformed_regions):
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

    def _check_overlap_batch(self, df, region_dict):
        omit_series = pd.Series(False, index=df.index)
        for chromosome, regions in region_dict.items():
            for start, end in regions:
                overlap_cond = (
                        ((df['chr_1'] == chromosome) & (df['start_1'] <= end) & (df['end_1'] >= start)) |
                        ((df['chr_2'] == chromosome) & (df['start_2'] <= end) & (df['end_2'] >= start))
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
