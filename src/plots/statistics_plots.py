# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from matplotlib import pyplot as plt
import dask.dataframe as dd

class PValuePlotter:

    def __init__(self, data: dd.DataFrame, config, metadata):
        self.data = data
        self.config = config
        self.metadata = metadata

    def plot_pvalues_vs_distance(self):
        sample_data = self.data.sample(frac=1).compute() if isinstance(self.data, dd.DataFrame) else self.data

        plt.figure(figsize=(10, 6))
        plt.scatter(sample_data["genomic_distance"], sample_data["p_value"], alpha=0.5)
        plt.title("P-Values vs Genomic Distance")
        plt.xlabel("Genomic Distance")
        plt.ylabel("P-Value")
        plt.xscale("log")
        plt.yscale("log")
        plt.grid(True, which="both", ls="--")
        plt.savefig(self.config.paths.plot_dir / f"{self.metadata.experiment}_{self.metadata.resolution}_pvalues_vs_distance_distribution.png")

    def plot_pvalue_distribution(self):
        sample_data = self.data.sample(frac=1).compute() if isinstance(self.data, dd.DataFrame) else self.data

        plt.figure(figsize=(10, 6))
        plt.hist(sample_data["p_value"], bins=100, alpha=0.5)
        plt.title("P-Value Distribution")
        plt.xlabel("P-Value")
        plt.ylabel("Frequency")
        plt.yscale("log")
        plt.grid(True, which="both", ls="--")
        plt.savefig(self.config.paths.plot_dir / f"{self.metadata.experiment}_{self.metadata.resolution}_pvalue_distribution.png")

    def plot_qvalue_distribution(self):
        sample_data = self.data.sample(frac=1).compute() if isinstance(self.data, dd.DataFrame) else self.data

        plt.figure(figsize=(10, 6))
        plt.hist(sample_data["q_value"], bins=100, alpha=0.5)
        plt.title("Q-Value Distribution")
        plt.xlabel("Q-Value")
        plt.ylabel("Frequency")
        plt.yscale("log")
        plt.grid(True, which="both", ls="--")
        plt.savefig(self.config.paths.plot_dir / f"{self.metadata.experiment}_{self.metadata.resolution}_qvalue_distribution.png")
