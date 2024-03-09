# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.containers import FilteringOutput, Metadata
from src.filtering.chromosome_filter import FilterChromosomes
from src.filtering.range_filter import FilterRanges
from src.filtering.filtering_utils import SplitByInteractionType, FilterBias
from src.filtering.blacklist_filter import RemoveBlacklistedRegions
from src.filtering.cytobands_filter import RemoveCytobandRegions
from src.data_preparation.interaction_calculator import InteractionCalculator
from src.filtering.filtering_utils import chromosome_key_sort
from typing import List
import dask.dataframe as dd

class FilteringController:

    def __init__(self, config: Config, interaction_output):
        self.config = config
        self.data: dd.DataFrame = interaction_output.data
        self.metadata = interaction_output.metadata
        self.blacklist_filter = RemoveBlacklistedRegions(config)
        self.cytoband_filter = RemoveCytobandRegions(config)

    def run(self) -> List[FilteringOutput]:
        # Split datasets by interaction type (inter and intra)
        splitter = SplitByInteractionType(self.data)
        intra_df, inter_df = splitter.split_datasets_by_interaction_type()

        # Instantiate metadata and output for intra and inter data
        intra_metadata = self._instantiate_metadata("intra")
        inter_metadata = self._instantiate_metadata("inter")
        intra_output = FilteringOutput(data=intra_df, metadata=intra_metadata)
        inter_output = FilteringOutput(data=inter_df, metadata=inter_metadata)

        # Filter outputs based on interaction_type specified in the config
        outputs = [intra_output, inter_output]
        if self.config.pipeline_settings.interaction_type != "mixed":
            outputs = [output for output in outputs if output.metadata.interaction_type == self.config.pipeline_settings.interaction_type]

        # Apply filters and return the filtered outputs
        for output in outputs:
            self._apply_filters(output)

        return outputs

    def _apply_filters(self, output):
        if self.config.statistical_settings.use_filtered_data_for_average_contact_probability:
            self._apply_chromosome_filtering(output)
            self._apply_range_filtering(output)
            self._apply_blacklist_filtering(output)
            self._apply_cytoband_filtering(output)
            self._apply_bias_filtering(output)
            self._calculate_interaction_counts(output)
        else:
            self._apply_chromosome_filtering(output)
            self._calculate_interaction_counts(output)
            self._apply_range_filtering(output)
            self._apply_blacklist_filtering(output)
            self._apply_cytoband_filtering(output)
            self._apply_bias_filtering(output)

    def _apply_chromosome_filtering(self, output: dd.DataFrame):
        if self.config.pipeline_settings.remove_chromosomes or self.config.pipeline_settings.select_chromosomes:
            chromosome_filter = FilterChromosomes(self.config)
            output.data = chromosome_filter.filter_data(output.data, output.metadata.interaction_type)

    def _apply_range_filtering(self, output: dd.DataFrame):
        range_filter = FilterRanges(self.config)
        output.data = range_filter.apply_filters(output.data)

    def _apply_blacklist_filtering(self, output: dd.DataFrame):
        if self.config.pipeline_settings.filter_blacklist:
            blacklist_filter = RemoveBlacklistedRegions(self.config)
            output.data = blacklist_filter.filter_blacklist(output.data, output.metadata.resolution)

    def _apply_cytoband_filtering(self, output: dd.DataFrame):
        if self.config.pipeline_settings.filter_cytobands:
            output.data = self.cytoband_filter.filter_cytobands(output.data, output.metadata.resolution)

    def _apply_bias_filtering(self, output: dd.DataFrame):
        if self.config.statistical_settings.use_hicpro_bias:
            bias_filter = FilterBias()
            output.data = bias_filter.filter_bias(output.data)


    @staticmethod
    def _calculate_interaction_counts(output):
        interaction_calculator = InteractionCalculator(output.data)

        if output.metadata.interaction_type == "intra":
            total_possible_intra, intra_per_chromosome_df = interaction_calculator.calculate_total_possible_bins_intra()

            # Convert DataFrame to dictionary with chromosomes as keys and counts as values
            chromosome_counts_dict = intra_per_chromosome_df.set_index("chr_1")["possible_intra"].to_dict()

            # Sort the chromosomes order and re-create dictionary in sorted order
            sorted_chromosome_keys = chromosome_key_sort(chromosome_counts_dict.keys())
            sorted_chromosome_counts_dict = {chrom: chromosome_counts_dict[chrom] for chrom in sorted_chromosome_keys}

            # Update metadata with sorted dictionary and total count
            output.metadata.max_possible_interacting_bins_intra = total_possible_intra
            output.metadata.max_possible_interacting_bins_per_chromosome_intra = sorted_chromosome_counts_dict

        elif output.metadata.interaction_type == "inter":
            total_possible_inter = interaction_calculator.calculate_total_possible_bins_inter()
            output.metadata.max_possible_interacting_bins_inter = total_possible_inter

    def _instantiate_metadata(self, interaction_type):
        # Instantiate metadata based on interaction type, including placeholder for interaction counts.
        return Metadata(
            experiment=self.metadata.experiment,
            resolution=self.metadata.resolution,
            interaction_type=interaction_type,
            bias_file_path=self.metadata.bias_file_path,
            max_possible_interacting_bins_intra=None,
            max_possible_interacting_bins_inter=None,
            max_possible_interacting_bins_per_chromosome_intra=None,
        )
