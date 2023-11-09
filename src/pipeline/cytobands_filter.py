# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import BlacklistOutput, CytobandOutput
from src.reference_handling.reference_parser import ReferenceFileParser as ReferenceParser
import pybedtools as pbt
import pandas as pd

class CytobandHelper:
    @staticmethod
    def get_cytoband_regions(config: Config) -> pbt.BedTool:
        reference_parser = ReferenceParser(config)
        cytoband_bedtool = reference_parser.parse_cytoband_file()
        # Filter for regions containing 'acen' (centromeric regions)
        acen_regions = cytoband_bedtool.filter(lambda x: x[4] == "acen").saveas()
        return acen_regions

class RemoveCytobandRegions:
    def __init__(self, config: Config, blacklist_output: BlacklistOutput):
        self.config = config
        self.blacklist_output = blacklist_output
        self.blacklist_ddf = blacklist_output.blacklist_ddf
        self.cytoband_regions = CytobandHelper.get_cytoband_regions(config)
        self.window_size = blacklist_output.metadata.resolution

    def _filter_single_partition(self, partition: pd.DataFrame) -> pd.DataFrame:
        pbt_partition = pbt.BedTool.from_dataframe(partition)
        # Here, we want to make sure we use the window parameter correctly
        filtered_partition = pbt_partition.window(self.cytoband_regions, w=self.window_size, v=True)
        filtered_partition_df = filtered_partition.to_dataframe()
        filtered_partition_df.columns = self.blacklist_output.blacklist_ddf.columns
        return filtered_partition_df

    def filter_cytobands(self):
        filtered_partitions = self.blacklist_ddf.map_partitions(
            self._filter_single_partition,
            meta=self.blacklist_ddf._meta  # use df struct from blacklist class (which again uses from bedpe class, ignore warning)
        )
        return filtered_partitions

    def run(self) -> CytobandOutput:
        # Triggers the actual computation
        filtered_ddf = self.filter_cytobands()
        result = filtered_ddf.compute()
        output = CytobandOutput(metadata=self.blacklist_output.metadata, cytoband_ddf=result)
        return output
