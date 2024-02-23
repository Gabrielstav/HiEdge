# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import FilteringOutput, Metadata
from src.filtering.chromosome_filter import FilterChromosomes
from src.filtering.range_filter import FilterRanges
from src.filtering.filtering_utils import SplitByInteractionType, FilterBias
from src.filtering.blacklist_filter import RemoveBlacklistedRegions
from src.filtering.cytobands_filter import RemoveCytobandRegions
from src.filtering.filtering_utils import chromosome_key_sort
from typing import List
import dask.dataframe as dd


class FilteringController:

    def __init__(self, config: Config, interaction_output):
        self.config = config
        self.data = interaction_output.data
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
        # Apply all the filters on the output
        self._apply_chromosome_filtering(output)
        self._apply_range_filtering(output)
        self._apply_blacklist_filtering(output)
        self._apply_cytoband_filtering(output)
        self._apply_bias_filtering(output)
        self._update_chromosomes_present_in_metadata(output)
        self._remove_metadata_for_filtered_chromosomes(output.metadata)


    def _split_datasets_by_interaction_type(self):
        splitter = SplitByInteractionType(self.data)
        self.intra_df, self.inter_df = splitter.split_datasets_by_interaction_type()
        self.intra_metadata = self._instantiate_metadata("intra")
        self.inter_metadata = self._instantiate_metadata("inter")

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

    def _instantiate_metadata(self, interaction_type):
        return Metadata(experiment=self.metadata.experiment,
                        resolution=self.metadata.resolution,
                        interaction_type=interaction_type,
                        bias_file_path=self.metadata.bias_file_path)

    def _remove_metadata_for_filtered_chromosomes(self, intra_metadata):
        if intra_metadata.interaction_count_per_chromosome_intra and self.config.statistical_settings.normalize_expected_frequency:
            intra_metadata.interaction_count_per_chromosome_intra = {
                k: v for k, v in intra_metadata.interaction_count_per_chromosome_intra.items()
                if k in intra_metadata.chromosomes_present
            }

    @staticmethod
    def _update_chromosomes_present_in_metadata(output):
        if output.metadata.interaction_type == "intra":
            unique_chromosomes_intra = output.data["chr_1"].unique().compute().tolist()
            chromosome_key_sort(unique_chromosomes_intra)
            output.metadata.chromosomes_present = unique_chromosomes_intra

        if output.metadata.interaction_type == "inter":
            chromosomes_1 = output.data["chr_1"].unique().compute().tolist()
            chromosomes_2 = output.data["chr_2"].unique().compute().tolist()
            unique_chromosomes_inter = list(set(chromosomes_1 + chromosomes_2))
            chromosome_key_sort(unique_chromosomes_inter)
            output.metadata.chromosomes_present = unique_chromosomes_inter

    @staticmethod
    def _instantiate_output(data, metadata):
        return FilteringOutput(data=data, metadata=metadata)
