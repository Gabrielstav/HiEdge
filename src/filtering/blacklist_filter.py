# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.reference_handling.reference_parser import ReferenceFileParser as ReferenceParser
import pandas as pd
import pybedtools as pbt
import dask.dataframe as dd
from typing import Union

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

    def filter_single_partition(self, partition: pd.DataFrame, resolution) -> pd.DataFrame:
        pbt_partition = pbt.BedTool.from_dataframe(partition)
        filtered_partition = pbt_partition.window(self.blacklisted_regions, r=False, v=True, w=resolution)
        filtered_partition_df = filtered_partition.to_dataframe(names=["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count", "idx_1", "idx_2", "bias_1", "bias_2"], engine="python")
        return filtered_partition_df

    def filter_blacklist(self, bedpe_data: Union[dd.DataFrame, pd.DataFrame], resolution) -> dd.DataFrame:
        if isinstance(bedpe_data, pd.DataFrame):
            # Convert pandas DataFrame to Dask DataFrame
            bedpe_ddf = dd.from_pandas(bedpe_data, npartitions=1)
        else:
            bedpe_ddf = bedpe_data

        filtered_partitions = bedpe_ddf.map_partitions(self.filter_single_partition, resolution)
        print(filtered_partitions.compute().head())

        return filtered_partitions

class RemoveBlacklistedRegionsDask:
    # TODO: Implement filtering logic in Dask, make series of blacklisted regions and filter on each partition using overlap or midpoints in region
    pass
