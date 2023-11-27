# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import BlacklistOutput, CytobandOutput, BedpeOutput
from src.reference_handling.reference_parser import ReferenceFileParser as ReferenceParser
import pybedtools as pbt
import pandas as pd

class CytobandInputResolver:
    def __init__(self, bedpe_output, blacklist_output):
        self.bedpe_output = bedpe_output
        self.blacklist_output = blacklist_output

    def resolve_input(self):
        if isinstance(self.blacklist_output, BlacklistOutput):
            return self.blacklist_output
        elif isinstance(self.bedpe_output, BedpeOutput):
            return self.bedpe_output
        else:
            raise ValueError("Invalid data class type provided")

class CytobandHelper:
    @staticmethod
    def get_cytoband_regions(config: Config) -> pbt.BedTool:
        reference_parser = ReferenceParser(config)
        cytoband_bedtool = reference_parser.parse_cytoband_file()
        # Filter for regions containing "acen" (centromeric regions)
        acen_regions = cytoband_bedtool.filter(lambda x: x[4] == "acen").saveas()
        return acen_regions

class RemoveCytobandRegions:
    def __init__(self, config: Config, data_input):
        self.config = config
        self.data_ddf = data_input.bedpe_ddf
        self.metadata = data_input.metadata
        self.cytoband_regions = CytobandHelper.get_cytoband_regions(config)
        self.window_size = self.metadata.resolution

    def _filter_single_partition(self, partition: pd.DataFrame) -> pd.DataFrame:
        pbt_partition = pbt.BedTool.from_dataframe(partition)
        filtered_partition = pbt_partition.window(self.cytoband_regions, w=self.window_size, v=True)
        filtered_partition_df = filtered_partition.to_dataframe()
        filtered_partition_df.columns = self.data_ddf.columns
        return filtered_partition_df

    def filter_cytobands(self):
        filtered_partitions = self.data_ddf.map_partitions(
            self._filter_single_partition,
            meta=self.data_ddf._meta
        )
        return filtered_partitions

    def run(self) -> CytobandOutput:
        filtered_ddf = self.filter_cytobands()
        result = filtered_ddf.compute()
        output = CytobandOutput(metadata=self.metadata, bedpe_ddf=result)
        return output
