# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config

class InteractionCalculator:
    def __init__(self, config: Config, data):
        self.config = config
        self.data = data  # interaction dask dataframe

    def calculate_intra_interactions(self):
        """
        Calculate total possible intrachromosomal interactions and interactions per chromosome for intra.
        """
        # chr_1 is used for intra calculations
        chrom_counts = self.data["chr_1"].value_counts().compute()

        total_bins = chrom_counts.sum()
        intra_counts = (chrom_counts * (chrom_counts + 1)) // 2
        total_possible_intra = intra_counts.sum()

        return total_possible_intra, intra_counts.to_dict()

    def calculate_inter_interactions(self):
        """
        Calculate total possible interchromosomal interactions and interactions per chromosome for inter.
        """
        # chr_1 and chr_2 are used for inter calculations
        chrom_counts_1 = self.data["chr_1"].value_counts().compute()
        chrom_counts_2 = self.data["chr_2"].value_counts().compute()
        combined_chrom_counts = chrom_counts_1.add(chrom_counts_2, fill_value=0)

        total_bins = combined_chrom_counts.sum()
        inter_counts = combined_chrom_counts * (total_bins - combined_chrom_counts)
        total_possible_inter = inter_counts.sum() // 2  # Divide by 2 to avoid double counting

        return total_possible_inter, combined_chrom_counts.to_dict()
