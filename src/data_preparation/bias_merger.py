# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
import dask.dataframe as dd

class BiasMerger:

    def __init__(self, config: Config, interaction_df, bias_file_path):
        self.config = config
        self.interaction_df = interaction_df
        self.bias_file_path = bias_file_path

    def process_and_merge_bias(self):
        bias_series = self._read_and_process_bias()
        return self.merge_bias(self.interaction_df, bias_series)

    def _read_and_process_bias(self):
        # read bias file into a Dask series
        bias_series = dd.read_csv(self.bias_file_path, header=None, squeeze=True)

        # replace "nan" string with -1 bias
        bias_series = bias_series.mask(bias_series == "nan", -1).astype(float)
        lower_bound = self.config.pipeline_settings.bias_lower_bound
        upper_bound = self.config.pipeline_settings.bias_upper_bound

        # filter out-of-bounds values
        bias_series = bias_series.where(bias_series.between(lower_bound, upper_bound), -1)
        return bias_series

    @staticmethod
    def merge_bias(merged_ddf, bias_series):
        # check if extending the bias series is necessary
        if len(bias_series) < (merged_ddf["id1"].max() + 1):

            # TODO: if mismatch in lengths are detected log which loci are missing
            # extend bias_series to match the length of bedpe_ddf (BED file length)
            extended_bias_series = bias_series.reindex(range(merged_ddf['id1'].max() + 1), fill_value=-1)
        else:
            # use original bias series if it covers all IDs
            extended_bias_series = bias_series

        # merge extended bias values with bedpe_ddf
        merged_ddf["bias_1"] = extended_bias_series.loc[merged_ddf["id1"]].values
        merged_ddf["bias_2"] = extended_bias_series.loc[merged_ddf["id2"]].values

        return merged_ddf
