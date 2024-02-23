# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from dask import dataframe as dd


class SplitByInteractionType:

    def __init__(self, bedpe_ddf: dd.DataFrame):
        self.bedpe_ddf = bedpe_ddf

    def split_datasets_by_interaction_type(self):
        intra_df = self.bedpe_ddf[self.bedpe_ddf["chr_1"] == self.bedpe_ddf["chr_2"]]
        inter_df = self.bedpe_ddf[self.bedpe_ddf["chr_1"] != self.bedpe_ddf["chr_2"]]
        return intra_df, inter_df


class FilterBias:
    @staticmethod
    def filter_bias(data: dd.DataFrame) -> dd.DataFrame:
        # Filter out interactions where either bias value is -1
        return data[(data["bias_1"] != -1) & (data["bias_2"] != -1)]


def chromosome_key_sort(chromosomes):
    custom_order = {"X": 24, "Y": 25, "M": 26}

    def sort_key(chromosome):
        num_part = "".join([char for char in chromosome if char.isdigit()])
        num_part = int(num_part) if num_part else 0

        # assign custom order to non-numeric chromosomes, default to 27 for others
        non_num_part = custom_order.get(chromosome.replace("chr", "").upper(), 27)

        return num_part, non_num_part

    sorted_chromosomes = sorted(chromosomes, key=sort_key)
    return sorted_chromosomes
