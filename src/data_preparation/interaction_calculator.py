# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config

class InteractionCalculator:
    def __init__(self, config: Config, data):
        self.config = config
        self.data = data  # interaction dask dataframe

    def calculate_total_interactions(self, data=None):
        chrom_counts = self.data["chr"].value_counts().compute()
        total_bins = chrom_counts.sum()

        intra_counts = (chrom_counts * (chrom_counts + 1)) // 2
        inter_counts = chrom_counts * (total_bins - chrom_counts)

        total_possible_intra = intra_counts.sum()
        total_possible_inter = inter_counts.sum() // 2

        return total_possible_intra, total_possible_inter, intra_counts.to_dict()
