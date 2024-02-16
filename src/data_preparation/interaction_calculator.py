# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config

class InteractionCalculator:
    def __init__(self, config: Config, data):
        self.config = config
        self.data = data  # interaction dask dataframe

    def calculate_total_interactions(self):
        chrom_counts_1 = self.data["chr_1"].value_counts().compute()
        chrom_counts_2 = self.data["chr_2"].value_counts().compute()

        # Combine counts from both chromosome columns, accounting for interactions involving different chromosomes
        chrom_counts = chrom_counts_1.add(chrom_counts_2, fill_value=0)

        total_bins = chrom_counts.sum()

        intra_counts = (chrom_counts * (chrom_counts + 1)) // 2
        inter_counts = chrom_counts * (total_bins - chrom_counts)

        total_possible_intra = intra_counts.sum()
        total_possible_inter = inter_counts.sum() // 2

        return total_possible_intra, total_possible_inter, chrom_counts.to_dict()
