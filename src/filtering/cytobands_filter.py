# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.reference_handling.reference_parser import ReferenceFileParser as ReferenceParser
import pybedtools as pbt
import pandas as pd
import dask.dataframe as dd

# mb use bedpe schema when returning data

class CytobandHelper:
    @staticmethod
    def get_cytoband_regions(config: Config) -> pbt.BedTool:
        reference_parser = ReferenceParser(config)
        cytoband_bedtool = reference_parser.parse_cytoband_file()
        # Filter for regions containing "acen" (centromeric regions)
        acen_regions = cytoband_bedtool.filter(lambda x: x[4] == "acen").saveas()
        return acen_regions

class RemoveCytobandRegions:

    def __init__(self, config: Config):
        self.config = config
        self.cytoband_regions = CytobandHelper.get_cytoband_regions(config)

    def filter_single_partition(self, partition: pd.DataFrame, resolution) -> pd.DataFrame:
        pbt_partition = pbt.BedTool.from_dataframe(partition)
        filtered_partition = pbt_partition.window(self.cytoband_regions, r=False, v=True, w=resolution)
        filtered_partition_df = filtered_partition.to_dataframe(names=["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count", "idx_1", "idx_2", "bias_1", "bias_2"])
        return filtered_partition_df

    def filter_cytobands(self, bedpe_ddf: dd.DataFrame, resolution) -> dd.DataFrame:
        filtered_partitions = bedpe_ddf.map_partitions(self.filter_single_partition, resolution)
        if not isinstance(filtered_partitions, dd.DataFrame):
            print("Filtered partitions are not a dask dataframe!")
        return filtered_partitions


class RemoveCytobandRegionsDask:
    # TODO: Implement filtering logic in Dask, make seris of cytoband regions and filter on each partition using overlap or midpoints in region
    pass
