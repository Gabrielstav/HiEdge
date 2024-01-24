# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import SplineInput
import pandas as pd
import dask.dataframe as dd

class FDRCalculator:

    def __init__(self, config, data):
        self.config = config
        self.data = data
        self.metadata = data.metadata
        self.p_value_column = "p_value"

    def run(self):

        # Intra data and partition method enabled
        if self.metadata.interaction_type == "intra" and self.config.use_partition_method_intra:
            return self._partition_based_fdr(self.config.pipeline_settings.num_partitions_intra)

        # Inter data and partition method enabled
        elif self.metadata.interaction_type == "inter" and self.config.use_partition_method_inter:
            return self._partition_based_fdr(self.config.pipeline_settings.num_partitions_inter)

        else:
            return self._sequential_fdr()

    def _partition_based_fdr(self, num_partitions):
        data_partitioned = self.data.map_partitions(
            lambda df: pd.qcut(df[self.p_value_column].rank(method="first"), num_partitions, labels=False),
            meta=("partition", "int64")
        )
        adjusted_partitions = data_partitioned.groupby("partition").apply(
            self._apply_bh_within_partition,
            meta=data_partitioned
        )
        return adjusted_partitions.drop("partition", axis=1).reset_index(drop=True)

    def _apply_bh_within_partition(self, df):
        df_sorted = df.sort_values(self.p_value_column)
        m = len(df_sorted)
        df_sorted["adjusted_p"] = [
            min(p * m / rank, 1) for rank, p in enumerate(df_sorted[self.p_value_column], start=1)
        ]
        return df_sorted

    def _sequential_fdr(self):
        collected_data = self.data.compute()
        sorted_data = collected_data.sort_values(self.p_value_column)
        m = len(sorted_data)
        sorted_data["adjusted_p"] = [
            min(p * m / rank, 1) for rank, p in enumerate(sorted_data[self.p_value_column], start=1)
        ]
        return sorted_data

    def _get_fdr_threshold(self):
        if self.metadata.interaction_type == "intra":
            return self.config.fdr_threshold_intra
        elif self.metadata.interaction_type == "inter":
            return self.config.fdr_threshold_inter
        else:
            raise ValueError("Unknown interaction type ")  # TODO: Remove this after making config validation

    @staticmethod
    def _apply_fdr_threshold(fdr_results, fdr_threshold):
        return fdr_results[fdr_results["adjusted_p"] <= fdr_threshold]


