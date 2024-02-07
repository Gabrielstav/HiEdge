# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import FilteringOutput, Metadata
from src.filtering.chromosome_filter import FilterChromosomes
from src.filtering.range_filter import FilterRanges
from src.filtering.filtering_utils import SplitByInteractionType, FilterInteractionDistances, FilterBias
from src.filtering.blacklist_filter import RemoveBlacklistedRegions
from src.filtering.cytobands_filter import RemoveCytobandRegions
from typing import List


class FilteringController:

    def __init__(self, config: Config, interaction_output):
        self.config = config
        # Assuming we get instance of dataclass containing data and metadata (duck typing)
        self.data = interaction_output.data
        self.metadata = interaction_output.metadata
        self.blacklist_filter = RemoveBlacklistedRegions(config)
        self.cytoband_filter = RemoveCytobandRegions(config)

    def run_filters(self) -> List[FilteringOutput]:
        # Split datasets by interaction type (inter and intra)
        splitter = SplitByInteractionType(self.data)
        intra_df, inter_df = splitter.split_datasets_by_interaction_type()

        # Instantiate metadata and output for intra and inter data
        intra_metadata = self._instantiate_metadata("intra")
        inter_metadata = self._instantiate_metadata("inter")
        intra_output = FilteringOutput(data=intra_df, metadata=intra_metadata)
        inter_output = FilteringOutput(data=inter_df, metadata=inter_metadata)

        # Apply the filters on intra_output and inter_output
        self._apply_chromosome_filtering(intra_output, inter_output)
        self._update_chromosomes_present_metadata(intra_metadata, inter_metadata)
        self._apply_range_filtering(intra_output, inter_output)
        self._apply_interaction_distance_filtering(intra_output, inter_output)
        self._apply_blacklist_filtering(intra_output, inter_output)
        self._apply_cytoband_filtering(intra_output, inter_output)
        self._apply_bias_filtering(intra_output, inter_output)

        return [intra_output, inter_output]

    def _split_datasets_by_interaction_type(self):
        splitter = SplitByInteractionType(self.data)
        self.intra_df, self.inter_df = splitter.split_datasets_by_interaction_type()
        self.intra_metadata = self._instantiate_metadata("intra")
        self.inter_metadata = self._instantiate_metadata("inter")

    def _apply_chromosome_filtering(self, intra_output, inter_output):
        if self.config.pipeline_settings.remove_chromosomes or self.config.pipeline_settings.select_chromosomes:
            chromosome_filter = FilterChromosomes(self.config)
            intra_output.data = chromosome_filter.filter_data(intra_output.data)
            inter_output.data = chromosome_filter.filter_data(inter_output.data)

    def _apply_range_filtering(self, intra_output, inter_output):
        range_filter = FilterRanges(self.config)
        if self.config.pipeline_settings.omit_regions:
            intra_output.data = range_filter.filter_omit_regions(intra_output.data)
            inter_output.data = range_filter.filter_omit_regions(inter_output.data)
        elif self.config.pipeline_settings.select_regions:
            intra_output.data = range_filter.filter_select_regions(intra_output.data)
            inter_output.data = range_filter.filter_select_regions(inter_output.data)

    def _apply_interaction_distance_filtering(self, intra_output, inter_output):
        if self.config.pipeline_settings.min_interaction_range or self.config.pipeline_settings.max_interaction_range:
            interaction_distance_filter = FilterInteractionDistances(self.config)
            intra_output.data = interaction_distance_filter.filter_data(intra_output.data)
            inter_output.data = interaction_distance_filter.filter_data(inter_output.data)

    def _apply_blacklist_filtering(self, intra_output, inter_output):
        if self.config.pipeline_settings.filter_blacklist:
            blacklist_filter = RemoveBlacklistedRegions(self.config)
            intra_output.data = blacklist_filter.filter_blacklist(intra_output.data)
            inter_output.data = blacklist_filter.filter_blacklist(inter_output.data)

    def _apply_cytoband_filtering(self, intra_output, inter_output):
        if self.config.pipeline_settings.filter_cytobands:
            intra_output.data = self.cytoband_filter.filter_cytobands(intra_output.data)
            inter_output.data = self.cytoband_filter.filter_cytobands(inter_output.data)

    def _apply_bias_filtering(self, intra_output, inter_output):
        if self.config.pipeline_settings.use_hicpro_bias:
            bias_filter = FilterBias()
            intra_output.data = bias_filter.filter_bias(intra_output.data)
            inter_output.data = bias_filter.filter_bias(inter_output.data)

    def _instantiate_metadata(self, interaction_type):
        return Metadata(experiment=self.metadata.experiment,
                        resolution=self.metadata.resolution,
                        interaction_type=interaction_type,
                        bias_file_path=self.metadata.bias_file_path)

    def _update_chromosomes_present_metadata(self, intra_metadata, inter_metadata):
        if self.config.pipeline_settings.remove_chromosomes:
            all_chromosomes = set(intra_metadata.interaction_count_per_chromosome.keys())
            removed_chromosomes = set(self.config.pipeline_settings.remove_chromosomes)
            remaining_chromosomes = list(all_chromosomes - removed_chromosomes)
            intra_metadata.chromosomes_present = remaining_chromosomes
            inter_metadata.chromosomes_present = remaining_chromosomes
        elif self.config.pipeline_settings.select_chromosomes:
            intra_metadata.chromosomes_present = self.config.pipeline_settings.select_chromosomes
            inter_metadata.chromosomes_present = self.config.pipeline_settings.select_chromosomes

    @staticmethod
    def _instantiate_output(data, metadata):
        return FilteringOutput(data=data, metadata=metadata)
