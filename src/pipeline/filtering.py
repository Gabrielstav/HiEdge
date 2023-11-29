# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.pipeline.bedpe_creator import BedpeOutput
from src.setup.data_structures import Metadata, SplittingOutput, FilteringOutput
from dataclasses import dataclass
from dask import dataframe as dd
import pandas as pd


# TODO: Make filtering on genomic ranges and chromosomes here

# TODO: Also make metadata updates on interaction type either at node level (?) or just on dataset level at least, we need this metadata for statistical testing
#   Either we split inter and intra here, before any processing, or we splint intra and inter before input to modeling..?
#   If we split before processing we could have intra and inter as separate instances and check for metadata attribute before passing to modeling
#   Or we could let the config decide, either process them together, intra only, inter only or process them separately...
#   If processed together, we could join inter and intra, reeturning unique entries but this could be kinda bad becasue different statistical assumptions in same dataset...
#   Since we store the dask dataframes containing interactions with metadata as fields in dataclasses, where one instance is unique data (data being one resolution + cell line),
#   We could update the dataclasses to also store interaction metadata (inter or intra) and then one instance could be interaction type.
#   The downside with this is that for one pipeline run, the intra data processing settings would also appply to inter data instances as the config is static for each run
#   Another appraoch would be to separate inter and intra settings in the config, and let users decide what resolution to process for both inter and intra, or just
#   have one option (inter and intra), always handling inter and intra in separate runs..
#   Perhaps the easiest thing, for the users, is to have one option in the config: process_intra, and then if this is true, set the resolutions at which the inter data should be processed.
#   For all resolutions not set by this option, only intra interactions would be considered. In the inter datasets, only inter interactions are considered. This way the dataclasses could just
#   have a interaction type field that is updated based on the interaction type in the data and the config settinggs, and since the instances are separate, we can just check the config setting to
#   only process inter data of certion resolutions and check the dataclasses when doing statistical stuff to know what method to apply to the data?

# TODO:
#   The simplest approach is to only let interaction type data be run-specific, it would be easier for config and statistical modeling
#   But, we could acieve the same thing with more complicated config, where inter and intra resolutions are set separately.
#   All intra datasets should then have inter interactions filtered out, and vice versa, where inter and intra resolutions are explicitly set in the config.
#   This has the advantage of not needing separate runs to process inter data, while still keeping the configuration flexible.
#   The downsie to splitting and letting the pipeline process both inter and intra is it makes the config so much more complicated.
#   When we create the bedpe daaframes, the containg both inter and intra interactions


class RunFiltering:

    def __init__(self):
        pass

    def on_data(self):
        """
        Run filtering classes on the datasets (bedpe dataclass containing the bedpe dask dataframes), instantiate dataclass "FilteringOutput" with data and metadata fields
        :return: Dataclass containing data and metadata for all unique sets of data with filtering applied
        """
        pass

class SplitInteractions:

    def __init__(self, config: Config, bedpe_output: BedpeOutput):
        self.config = config
        self.bedpe_metadata = bedpe_output.metadata
        self.bedpe_ddf = bedpe_output.bedpe_ddf

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

    def _check_overlap(self, interaction, region_dict):
        def overlaps_with(locus_chr, locus_start, locus_end):
            for start, end in region_dict.get(locus_chr, []):
                if locus_start <= end and locus_end >= start:
                    return True
            return False

        if interaction['interaction_type'] == 'intra':
            return overlaps_with(interaction['chr_1'], interaction['start_1'], interaction['end_1'])
        elif interaction['interaction_type'] == 'inter':
            return overlaps_with(interaction['chr_1'], interaction['start_1'], interaction['end_1']) or \
                   overlaps_with(interaction['chr_2'], interaction['start_2'], interaction['end_2'])
        return False

    def filter_omit_regions(self, data: dd.DataFrame) -> dd.DataFrame:
        omit_dict = self._get_omit_dict()
        filtered_data = data.map_partitions(
            lambda df: df[~df.apply(lambda row: self._check_overlap(row, omit_dict), axis=1)]
        )
        return filtered_data

    def filter_select_regions(self, data: dd.DataFrame) -> dd.DataFrame:
        select_dict = self._get_select_dict()
        filtered_data = data.map_partitions(
            lambda df: df[df.apply(lambda row: self._check_overlap(row, select_dict), axis=1)]
        )
        return filtered_data

# TODO: Finish Filtering on chromosomes

class FilterChromosomes:

    def __init__(self, config: Config):
        self.config = config

    def remove_chromosomes(self):
        """
        Removes specifie chromosomes from dataset
        :return: Dask dataframe in bedpe format not containing removed chromosomes
        """
        pass

    def include_chromosomes(self):
        """
        Removes all chromosomes not specified to be inlcuded
        :return: Dask dataframe in bedpe format only containing chromosomes included
        """
        pass