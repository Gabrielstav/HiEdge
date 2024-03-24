# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from matplotlib import pyplot as plt
import dask.dataframe as dd

# TODO: Improve this later, add more descriptive plots, log scale, increase resolution, implement the other methods.

class SplinePlotter:
    def __init__(self, spline_fitter, metadata, config):
        self.config = config
        self.metadata = metadata
        self.spline_fitter = spline_fitter

    def plot_spline_fit(self):
        x = self.spline_fitter.x
        y = self.spline_fitter.y
        y_mono = self.spline_fitter.mono_spline(x)

        plt.figure(figsize=(10, 6))
        plt.scatter(x, y, alpha=0.5, label="Data")
        plt.plot(x, y_mono, color="red", label="Fitted Spline")
        plt.title(f"Genomic Distance vs Probability in {self.metadata.experiment} {self.metadata.resolution}")
        plt.xlabel("Genomic Distance")
        plt.ylabel("Probability")
        plt.legend()
        plt.savefig(self.config.paths.plot_dir / f"{self.metadata.experiment}_{self.metadata.resolution}_spline_fit.png")

    def plot_spline_fit_with_error(self):
        pass

    def plot_spline_fit_with_confidence_intervals(self):
        pass

    def plot_distance_distributuion(self):
        x = self.spline_fitter.x
        plt.figure(figsize=(10, 6))
        plt.scatter(x, range(len(x)), alpha=0.5)
        plt.title(f"Distance Distribution in {self.metadata.experiment} {self.metadata.resolution}")
        plt.xlabel("Genomic Distance")
        plt.ylabel("Frequency")
        plt.savefig(self.config.paths.plot_dir / f"{self.metadata.experiment}_{self.metadata.resolution}_distance_distribution.png")

    def plot_probability_distribution(self):
        y = self.spline_fitter.y
        plt.figure(figsize=(10, 6))
        plt.scatter(y, range(len(y)), alpha=0.5)
        plt.title(f"Probability Distribution of metabins in {self.metadata.experiment} {self.metadata.resolution}")
        plt.xlabel("Probability")
        plt.ylabel("Frequency")
        plt.savefig(self.config.paths.plot_dir / f"{self.metadata.experiment}_{self.metadata.resolution}_probability_distribution.png")


class DistanceDecayPlotter:

    def __init__(self, interaction_data: dd.DataFrame, metadata, config):
        self.config = config
        self.metadata = metadata
        self.interaction_data = interaction_data

    def plot_distance_dependant_decay(self):
        x = self.interaction_data["genomic_distance"].values
        y = self.interaction_data["interaction_count"].values

        plt.figure(figsize=(10, 6))
        plt.title(f"Distance Dependant Decay in {self.metadata.experiment} {self.metadata.resolution}")
        plt.xlabel("Genomic Distance")
        plt.ylabel("Contact Count")
        plt.scatter(x, y, alpha=0.5)
        plt.savefig(self.config.paths.plot_dir / f"{self.metadata.experiment}_{self.metadata.resolution}_distance_decay.png")
