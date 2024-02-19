# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config

class InteractionCalculator:
    def __init__(self, data):
        self.data = data  # interaction dask dataframe

    def calculate_intra_interactions(self):
        chrom_counts = self.data["chr_1"].value_counts().compute()
        total_possible_intra = 0
        intra_counts_per_chromosome = {}

        for chrom, count in chrom_counts.items():
            intra_counts = (count * (count + 1)) // 2
            total_possible_intra += intra_counts
            intra_counts_per_chromosome[chrom] = intra_counts

        return total_possible_intra, intra_counts_per_chromosome

    def calculate_inter_interactions(self):
        chrom_counts = {chrom: self.data[self.data["chr_1"] == chrom]["chr_1"].count().compute() for chrom in self.data["chr_1"].unique()}
        total_bins = sum(chrom_counts.values())
        total_possible_inter = 0

        for chrom, n in chrom_counts.items():
            # Subtract current chromosome's bin count from total to calculate interchromosomal pairs
            inter_count_for_chrom = n * (total_bins - n)
            total_possible_inter += inter_count_for_chrom

        # Divide by 2 to avoid double counting
        total_possible_inter /= 2

        return total_possible_inter
