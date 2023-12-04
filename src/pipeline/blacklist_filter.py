# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import BedpeOutput, BlacklistOutput
from src.reference_handling.reference_parser import ReferenceFileParser as ReferenceParser
import pandas as pd
import pybedtools as pbt

class BlacklistHelper:
    @staticmethod
    def get_blacklisted_regions(config: Config) -> pbt.BedTool:
        reference_parser = ReferenceParser(config)
        return reference_parser.parse_blacklist_file()

class RemoveBlacklistedRegions:
    def __init__(self, config: Config, bedpe_output: BedpeOutput):
        self.config = config
        self.bedpe_output = bedpe_output
        self.bedpe_ddf = bedpe_output.data
        self.blacklisted_regions = BlacklistHelper.get_blacklisted_regions(config)
        self.window_size = bedpe_output.metadata.resolution

    def _filter_single_partition(self, partition: pd.DataFrame) -> pd.DataFrame:
        pbt_partition = pbt.BedTool.from_dataframe(partition)
        filtered_partition = pbt_partition.window(self.blacklisted_regions, v=True, wa=True)
        filtered_partition_df = filtered_partition.to_dataframe()
        filtered_partition_df.columns = self.bedpe_output.data.columns
        return filtered_partition_df

    def filter_blacklist(self):
        filtered_partitions = self.bedpe_ddf.map_partitions(
            self._filter_single_partition,
            meta=self.bedpe_ddf._meta # use df struct from bedpe class (ignore warning, Dask uses private naming convention for non-private methods)
        )
        return filtered_partitions

    def run(self) -> BlacklistOutput:
        filtered_ddf = self.filter_blacklist()
        result = filtered_ddf.compute()
        output = BlacklistOutput(metadata=self.bedpe_output.metadata, bedpe_ddf=result)
        return output
