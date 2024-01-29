# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.reference_handling.reference_parser import ReferenceFileParser as ReferenceParser
import pandas as pd
import pybedtools as pbt
import dask.dataframe as dd

# mb use bedpe schema when returning data

class BlacklistHelper:
    @staticmethod
    def get_blacklisted_regions(config: Config) -> pbt.BedTool:
        reference_parser = ReferenceParser(config)
        return reference_parser.parse_blacklist_file()

class RemoveBlacklistedRegions:
    def __init__(self, config: Config):
        self.config = config
        self.blacklisted_regions = BlacklistHelper.get_blacklisted_regions(config)

    def filter_single_partition(self, partition: pd.DataFrame) -> pd.DataFrame:
        pbt_partition = pbt.BedTool.from_dataframe(partition)
        filtered_partition = pbt_partition.window(self.blacklisted_regions, v=True, wa=True)
        filtered_partition_df = filtered_partition.to_dataframe()
        return filtered_partition_df

    def filter_blacklist(self, bedpe_ddf: dd.DataFrame) -> dd.DataFrame:
        filtered_partitions = bedpe_ddf.map_partitions(self.filter_single_partition)
        return filtered_partitions
