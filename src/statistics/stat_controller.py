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
from src.plots.spline_plots import SplinePlotter, DistanceDecayPlotter
from src.plots.statistics_plots import PValuePlotter


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

        if self.config.pipeline_settings.make_plots:
            DistanceDecayPlotter(filtered_data, metadata=self.metadata, config=self.config).plot_distance_dependant_decay()

        metabinned_data_with_stats = IntraStatsProcessor(self.config, filtered_data, self.metadata.unique_bins_per_chromosome, self.metadata.resolution).run()

        spline = self._fit_spline(metabinned_data_with_stats)
        if self.config.pipeline_settings.make_plots:
            SplinePlotter(spline, metadata=self.metadata, config=self.config).plot_spline_fit()
            SplinePlotter(spline, metadata=self.metadata, config=self.config).plot_distance_distributuion()
            SplinePlotter(spline, metadata=self.metadata, config=self.config).plot_probability_distribution()

        pvals = self._calculate_intra_pvals(spline, filtered_data)
        if self.config.pipeline_settings.make_plots:
            PValuePlotter(pvals, metadata=self.metadata, config=self.config).plot_pvalues_vs_distance()
            PValuePlotter(pvals, metadata=self.metadata, config=self.config).plot_pvalue_distribution()

        fdr = self._calculate_fdr_intra(pvals)
        if self.config.pipeline_settings.make_plots:
            PValuePlotter(fdr, metadata=self.metadata, config=self.config).plot_qvalue_distribution()
        return StatisticalOutput(metadata=self.metadata, data=fdr)

    def _process_inter_stats(self):
        inter_stats = self._prepare_inter_data(self.data)
        pvals = self._process_inter_pvals(inter_stats)
        fdr = self._calculate_fdr_inter(pvals)
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

    @staticmethod
    def _fit_spline(metabinned_data):
        spline_fitter = SplineFitter(metabinned_data)
        spline_fitter.fit_spline()
        return spline_fitter

    # def _spline_stats(self, spline, spline_input):
    #     spline_analysis = SplineAnalysis(self.config, spline, spline_input)
    #     spline_stats = spline_analysis.calculate_mse(), spline_analysis.calculate_residuals(), spline_analysis.calculate_r_squared()
    #     return spline_stats

    def _calculate_intra_pvals(self, spline, filtered_data):
        pval_calculator = IntraPValueCalculator(filtered_data, spline, self.metadata, self.config, self.metadata.unique_bins_per_chromosome)
        pvals = pval_calculator.run()
        return pvals

    def _process_inter_pvals(self, pval_input):
        pval_calculator = InterPValueCalculator(pval_input, self.config, self.metadata, self.metadata.max_possible_interaction_count_inter)
        pvals = pval_calculator.run()
        return pvals

    def _calculate_fdr_intra(self, pvals):
        fdr_calculator = FDRCalculator(pvals, self.config, self.metadata)
        fdr = fdr_calculator.run()
        return fdr

    def _calculate_fdr_inter(self, pvals):
        fdr_calculator = FDRCalculator(pvals, self.config)
        fdr = fdr_calculator.run()
        return fdr
