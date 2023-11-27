# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import ConfigMapper as Config
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

# TODO: Finish Filtering on ranges (inter, parsing) and same for chromosome (ez) 

class FilterRanges:

    def __init__(self, config: Config, data_input: SplittingOutput):
        self.config = config
        self.data_df = data_input.data
        self.metadata = data_input.metadata

    def filter_data(self):
        """
        Filter out regions based on config file settings
        :return: Dataclass with data and metadata, where the dd.DataFrames are filtered based
        """
        if self.metadata.interaction_type == 'intra':
            self.data_df = self._filter_intra_ranges(self.data_df)
        elif self.metadata.interaction_type == 'inter':
            self.data_df = self._filter_inter_ranges(self.data_df)

        # Return updated data
        return FilteringOutput(data=self.data_df, metadata=self.metadata)

    def _parse_ranges(self):
        pass


    def _filter_intra_ranges(self, df):
        if self.config.select_regions:
            df = self.include_selected_regions(df, self.config.select_regions)
        if self.config.omit_regions:
            df = self.exclude_omitted_regions(df, self.config.omit_regions)
        return df

    @staticmethod
    def include_selected_regions(df, select_regions):
        conditions = []
        for chromosome, regions in select_regions.items():
            for region in regions:
                condition = ((df['chr_1'] == chromosome) &
                             (df['start_1'] >= region.start) & (df['end_1'] <= region.end) &
                             (df['chr_2'] == chromosome) &
                             (df['start_2'] >= region.start) & (df['end_2'] <= region.end))
                conditions.append(condition)

        combined_condition = pd.concat(conditions, axis=1).any(axis=1)
        return df[combined_condition]

    @staticmethod
    def exclude_omitted_regions(df, omit_regions):
        conditions = []
        for chromosome, regions in omit_regions.items():
            for region in regions:
                condition = ~((df['chr_1'] == chromosome) &
                              ((df['start_1'] <= region.end) & (df['end_1'] >= region.start)) |
                              (df['chr_2'] == chromosome) &
                              ((df['start_2'] <= region.end) & (df['end_2'] >= region.start)))
                conditions.append(condition)

        combined_condition = pd.concat(conditions, axis=1).all(axis=1)
        return df[combined_condition]

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