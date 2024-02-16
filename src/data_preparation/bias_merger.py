# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
import dask.dataframe as dd

class BiasMerger:

    def __init__(self, config: Config, interaction_df, bias_file_path, bed_length=None):
        self.config = config
        self.interaction_df = interaction_df
        self.bias_file_path = bias_file_path
        self.bed_length = bed_length

    def process_and_merge_bias(self):
        bias_series = self._read_and_process_bias()
        print(f"Before merging: {bias_series.head(5)}")

        # Convert bias_series into a DataFrame and reset the index to get 'idx' as a column
        bias_df = bias_series.to_frame(name="bias_value").reset_index().rename(columns={"index": "idx"})
        bias_df["idx"] += 1  # Set start index to match index of interaction DataFrame (from matrix id)

        # Merge bias values based on idx with interaction DataFrame
        self.interaction_df = self.merge_bias_with_interaction(self.interaction_df, bias_df, self.bed_length)
        print(f"After merging: {self.interaction_df.head(10)}")

        return self.interaction_df

    def _read_and_process_bias(self):
        # Read bias file into a Dask DataFrame
        bias_df = dd.read_csv(self.bias_file_path, sep="\t", header=None, names=["bias_value"], encoding="latin-1")

        # Replace "nan" string with -1 bias
        bias_df["bias_value"] = bias_df["bias_value"].mask(bias_df["bias_value"] == "nan", -1).astype(float)

        lower_bound = self.config.statistical_settings.bias_lower_bound
        upper_bound = self.config.statistical_settings.bias_upper_bound

        bias_df["bias_value"] = bias_df["bias_value"].where(bias_df["bias_value"].between(lower_bound, upper_bound, inclusive="both"), -1)

        bias_series = bias_df["bias_value"]
        return bias_series

    @staticmethod
    def merge_bias_with_interaction(interaction_ddf, bias_df, bed_length):
        # Ensure bias_df is extended to match bed_length if necessary
        if bed_length and len(bias_df) < bed_length:
            # Assuming 'idx' is a continuous range from 0 to bed_length-1
            extended_bias_df = bias_df.reindex(range(bed_length), fill_value=-1)
        else:
            extended_bias_df = bias_df

        # Apply bias values to the interaction DataFrame by merging on "idx"
        interaction_ddf = interaction_ddf.merge(extended_bias_df, left_on="idx_1", right_on="idx", how="left").rename(columns={"bias_value": "bias_1"})
        interaction_ddf = interaction_ddf.merge(extended_bias_df, left_on="idx_2", right_on="idx", how="left").rename(columns={"bias_value": "bias_2"})

        if "idx_x" in interaction_ddf.columns:
            interaction_ddf = interaction_ddf.drop(columns=["idx_x"])
        if "idx_y" in interaction_ddf.columns:
            interaction_ddf = interaction_ddf.drop(columns=["idx_y"])

        print(interaction_ddf.head(5))
        return interaction_ddf
