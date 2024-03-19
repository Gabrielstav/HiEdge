# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.containers import StatisticalOutput
from src.statistics.intra_statistics import IntraStatsProcessor
from src.statistics.stat_utils import MidpointCalculator, DistanceCalculator, DistanceFilter
from src.statistics.inter_statistics import InterStatisticsCalculator
from src.statistics.spline_fit import SplineFitter, SplineAnalysis
from src.statistics.pvals import IntraPValueCalculator, InterPValueCalculator
from src.statistics.fdr import FDRCalculator
from src.setup.config_loader import Config
from matplotlib import pyplot as plt


# TODO: After q-vals, determine outliers and re-run spline fitting based on num-pass from config with outliers removed (basically just re-run this moethod with outliers removed)

class StatController:

    def __init__(self, config: Config, filtering_output) -> None:
        self.config = config
        self.metadata = filtering_output.metadata
        self.data = filtering_output.data

    def run(self):
        if self.metadata.interaction_type == "intra":
            print(f"Running intra stats for experiment: {str(self.metadata.experiment)} with resolution: {str(self.metadata.resolution)}")
            print(f"Metadata in stat controller: {self.metadata}")
            return self._process_intra_stats()
        elif self.metadata.interaction_type == "inter":
            return self._process_inter_stats()
        else:
            raise ValueError("Invalid interaction type")

    def _process_intra_stats(self):
        filtered_data = self._filter_on_distances()

        # plot data before metabinning to check for distance dependant decay
        distance_before_metabinning = filtered_data["genomic_distance"]
        contact_probability_before_metabinning = filtered_data["interaction_count"]
        plot_data(distance_before_metabinning, contact_probability_before_metabinning, title='Genomic Distance vs. Contact Count', xlabel='Genomic Distance', ylabel='Contact Count')

        metabinned_data_with_stats = IntraStatsProcessor(self.config, filtered_data, self.metadata.unique_bins_per_chromosome, self.metadata.resolution).run()

        # plot data after metabinning to check for distance dependant decay
        distance_after_metabinning = metabinned_data_with_stats["average_genomic_distance"]
        contact_probability_after_metabinning = metabinned_data_with_stats["average_contact_probability"]
        plot_data(distance_after_metabinning, contact_probability_after_metabinning, title='Average Genomic Distance vs. Contact Probability', xlabel='Average Genomic Distance', ylabel='Average Contact Probability')

        spline = self._fit_spline(metabinned_data_with_stats)
        pvals = self._calculate_intra_pvals(spline, filtered_data)
        fdr = self._calculate_fdr_intra(pvals)
        return StatisticalOutput(metadata=self.metadata, data=fdr)

    def _process_inter_stats(self):
        inter_stats = self._prepare_inter_data(self.data)
        pvals = self._process_pvals_inter(inter_stats)
        fdr = self._process_fdr_inter(pvals)
        return StatisticalOutput(metadata=self.metadata, data=fdr)

    def _filter_on_distances(self):
        midpoints = MidpointCalculator(self.config).calculate_midpoints(self.data)
        distances = DistanceCalculator(self.config).calculate_distances(midpoints)
        if self.config.pipeline_settings.interaction_distance_filters and self.metadata.resolution in self.config.pipeline_settings.interaction_distance_filters:
            filtered_data = DistanceFilter(self.config, distances, resolution=self.metadata.resolution).filter_distances()
            return_data = filtered_data
        else:
            return_data = distances
        return return_data

    def _prepare_inter_data(self, inter_data):
        inter_stats = InterStatisticsCalculator(inter_data, self.config).compute_inter_stats()
        return inter_stats

    def _fit_spline(self, metabinned_data):
        spline_fitter = SplineFitter(metabinned_data, self.config)
        spline_fitter.fit_spline().enforce_monotonicity()
        # plot spline
        spline_fitter.plot_spline()
        return spline_fitter

    def _spline_stats(self, spline, spline_input):
        spline_analysis = SplineAnalysis(self.config, spline, spline_input)
        spline_stats = spline_analysis.calculate_mse(), spline_analysis.calculate_residuals(), spline_analysis.calculate_r_squared()
        return spline_stats

    def _calculate_intra_pvals(self, spline, filtered_data):
        pval_calculator = IntraPValueCalculator(filtered_data, spline, self.metadata, self.config)
        pvals = pval_calculator.run()
        return pvals

    def _process_pvals_inter(self, pval_input):
        pval_calculator = InterPValueCalculator(pval_input, self.config, self.metadata, self.metadata.max_possible_interaction_count_inter)
        pvals = pval_calculator.run()
        return pvals

    def _calculate_fdr_intra(self, pvals):
        fdr_calculator = FDRCalculator(pvals, self.config)
        fdr = fdr_calculator.run()
        return fdr

    def _process_fdr_inter(self, pvals):
        fdr_calculator = FDRCalculator(pvals, self.config)
        fdr = fdr_calculator.run()
        return fdr

def plot_data(x, y, title='Average Genomic Distance vs. Contact Probability', xlabel='Average Genomic Distance', ylabel='Average Contact Probability'):
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, alpha=0.5)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
