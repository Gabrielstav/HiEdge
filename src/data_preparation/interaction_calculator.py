# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import dask.dataframe as dd

class InteractionCalculator:
    def __init__(self, data: dd.DataFrame):
        self.data: dd.DataFrame = data

        if not isinstance(self.data, dd.DataFrame):
            print(f"Converting to dask dataframe in InteractionCalculator: {type(self.data)}")
            self.data = dd.from_pandas(self.data, npartitions=1)

    def calculate_total_interactions_intra(self):
        # sum of interaction counts for all intra-chromosomal interactions
        total_interactions_intra = self.data[self.data["chr_1"] == self.data["chr_2"]]["interaction_count"].sum().compute()
        return total_interactions_intra

    def calculate_total_interactions_inter(self):
        # sum of interaction counts for all inter-chromosomal interactions
        total_interactions_inter = self.data[self.data["chr_1"] != self.data["chr_2"]]["interaction_count"].sum().compute()
        return total_interactions_inter

    def calculate_total_interactions_per_chromosome_intra(self):
        # sum of interaction counts per chromosome for all intra-chromosomal interactions
        total_interactions_per_chromosome_intra = self.data[self.data["chr_1"] == self.data["chr_2"]].groupby("chr_1")["interaction_count"].sum().compute()
        return total_interactions_per_chromosome_intra

    def calculate_unique_possible_distances_intra_per_chromosome(self):
        # Calculate unique distances per chromosome
        unique_distances_per_chromosome = self.data[self.data["chr_1"] == self.data["chr_2"]].groupby("chr_1")["genomic_distance"].nunique().compute()
        print(f"Unique distances per chromosome: {unique_distances_per_chromosome}")
        return unique_distances_per_chromosome

class PossiblePairsCalculator:
    def __init__(self, data: dd.DataFrame):
        self.data: dd.DataFrame = data

        if not isinstance(self.data, dd.DataFrame):
            print(f"Converting to dask dataframe in InteractionCalculator: {type(self.data)}")
            self.data = dd.from_pandas(self.data, npartitions=1)

    def calculate_total_possible_bins_intra(self):
        # Count unique bins per chromosome
        unique_bins_per_chr = self.data[self.data["chr_1"] == self.data["chr_2"]][["chr_1", "start_1"]].drop_duplicates().groupby("chr_1").size().compute()

        # Calculate possible intra interactions per chromosome
        intra_per_chromosome = (unique_bins_per_chr * (unique_bins_per_chr + 1)) // 2

        # Sum up to get total possible intra interactions
        total_possible_intra = intra_per_chromosome.sum()

        # Create DataFrame for per chromosome interaction counts
        intra_per_chromosome_df = intra_per_chromosome.reset_index(name="possible_intra")

        return total_possible_intra, intra_per_chromosome_df

    def calculate_total_possible_bins_inter(self):
        # Count unique bins across the entire dataset
        total_unique_bins = self.data[["chr_1", "start_1"]].drop_duplicates().shape[0].compute()

        # Calculate total possible inter interactions
        total_possible_inter = (total_unique_bins * (total_unique_bins - 1)) // 2

        return total_possible_inter

    def calculate_possible_pairs_within_metabin_distance(self, metabin_dataframe: dd.DataFrame):
        # is this possible pairs within each metabin? yes. So similar to possible pairs within each chromosome but for each metabin, same calculation?
        # count unique bins in each metabin, then use these bins to calculate possible pairs.
        # since my data is already filtered for in-range interactions, isn't every pair by definition inside the metabin distance range?
        # or is this the proportion of the total possible pairs that are contained in one metabin?
        pass