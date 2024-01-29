# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import BlacklistOutput, CytobandOutput, FilteringOutput, BedpeOutput, Metadata
from src.filtering.filtering import SplitByInteractionType, FilterChromosomes, FilterRanges, FilterInteractionDistances, FilterBias
from src.filtering.blacklist_filter import RemoveBlacklistedRegions as FilterBlacklist
from src.filtering.cytobands_filter import RemoveCytobandRegions as FilterCytobands
from dask import dataframe as dd

# TODO:
#   Make the Blacklist and Cytoband filtering classes work with the FilteringController
#   The goal is that the controller always returns instances of FilteringOutput if any filtering is enabled in config, if not, it returns the BedpeOutput instances

class FilteringController:

    def __init__(self, config: Config, data: dd.DataFrame, metadata: Metadata):
        self.config = config
        self.data = data
        self.metadata = metadata
        self.blacklist_filter = FilterBlacklist(self.config)
        self.cytoband_filter = FilterCytobands(self.config)

    def run_filters(self):

        # Slit datasets by interaction type (inter and intra)
        splitter = SplitByInteractionType(self.data)
        intra_df, inter_df = splitter.split_datasets_by_interaction_type()

        # Instantiate metadata and output for intra and inter data
        intra_metadata = self._instantiate_metadata("intra")
        inter_metadata = self._instantiate_metadata("inter")
        intra_output = self._instantiate_output(intra_df, intra_metadata)
        inter_output = self._instantiate_output(inter_df, inter_metadata)

        # Apply chromosome filtering if specified
        if self.config.pipeline_settings.remove_chromosomes or self.config.pipeline_settings.select_chromosomes:
            chromosome_filter = FilterChromosomes(self.config)
            intra_output.data = chromosome_filter.filter_data(intra_output.data)
            inter_output.data = chromosome_filter.filter_data(inter_output.data)

        # Apply range filtering if specified
        range_filter = FilterRanges(self.config)
        if self.config.pipeline_settings.omit_regions:
            intra_output.data = range_filter.filter_omit_regions(intra_output.data)
            inter_output.data = range_filter.filter_omit_regions(inter_output.data)
        elif self.config.pipeline_settings.select_regions:
            intra_output.data = range_filter.filter_select_regions(intra_output.data)
            inter_output.data = range_filter.filter_select_regions(inter_output.data)

        # Apply interaction distance filtering if specified
        if self.config.pipeline_settings.min_interaction_range or self.config.pipeline_settings.max_interaction_range:
            interaction_distance_filter = FilterInteractionDistances(self.config)
            intra_output.data = interaction_distance_filter.filter_data(intra_output.data)
            inter_output.data = interaction_distance_filter.filter_data(inter_output.data)

        # Apply blacklist filtering if specified
        if self.config.pipeline_settings.filter_blacklist:
            intra_output.data = self.blacklist_filter.filter_blacklist(intra_output.data)
            inter_output.data = self.blacklist_filter.filter_blacklist(inter_output.data)

        # Apply cytoband filtering if specified
        if self.config.pipeline_settings.filter_cytobands:
            intra_output.data = self.cytoband_filter.filter_cytobands(intra_output.data)
            inter_output.data = self.cytoband_filter.filter_cytobands(inter_output.data)

        # Apply bias filtering if specified
        if self.config.pipeline_settings.use_hicpro_bias:
            bias_filter = FilterBias()
            intra_output.data = bias_filter.filter_bias(intra_output.data)
            inter_output.data = bias_filter.filter_bias(inter_output.data)

        # Return the potentially modified FilteringOutput instances
        return intra_output, inter_output

    def _instantiate_metadata(self, interaction_type):
        return Metadata(experiment=self.metadata.experiment,
                        resolution=self.metadata.resolution,
                        interaction_type=interaction_type,
                        bias_file_path=self.metadata.bias_file_path)

    @staticmethod
    def _instantiate_output(data, metadata):
        return FilteringOutput(data=data, metadata=metadata)
