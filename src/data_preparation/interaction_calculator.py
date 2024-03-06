# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import dask.dataframe as dd

class InteractionCalculator:
    def __init__(self, data: dd.DataFrame):
        self.data: dd.DataFrame = data

        if not isinstance(self.data, dd.DataFrame):
            # convert to dask dataframe
            self.data = dd.from_pandas(self.data, npartitions=1)

    def calculate_intra_interactions(self):
        # Calculate unique bins per chromosome
        unique_bins_per_chr = self.data.groupby("chr_1").apply(lambda df: df[["start_1"]].drop_duplicates().shape[0], meta="int").compute()

        # Calculate possible intra interactions per chromosome
        intra_per_chromosome = unique_bins_per_chr * (unique_bins_per_chr - 1) // 2

        # Sum up to get total possible intra interactions
        total_possible_intra = intra_per_chromosome.sum()

        # Create DataFrame for per chromosome interaction counts
        intra_per_chromosome_df = intra_per_chromosome.reset_index(name="possible_intra")

        print(f"Total possible intra interactions: {total_possible_intra}"
              f"\nIntra interactions per chromosome:\n{intra_per_chromosome_df}")

        return total_possible_intra, intra_per_chromosome_df

    def calculate_inter_interactions(self):
        # Count unique bins across the entire dataset
        total_unique_bins = self.data[["chr_1", "start_1"]].drop_duplicates().shape[0].compute()

        # Calculate total possible inter interactions
        total_possible_inter = total_unique_bins * (total_unique_bins + 1) // 2   # - self.calculate_intra_interactions()[0] no need for this if calling after splitting in filtering_controller

        return total_possible_inter

    def calculate_per_chromosome_interactions(self):
        # This method is relevant for intra-chromosomal interactions only
        # Similar to calculate_intra_interactions but keeps the DataFrame for per chromosome calculations
        _, intra_per_chromosome_df = self.calculate_intra_interactions()

        return intra_per_chromosome_df
