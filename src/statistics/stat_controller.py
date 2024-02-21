# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.data_structures import Metadata, StatisticalOutput
from src.statistics.intra_statistics import DistanceCalculator, EqualOccupancyBinner, FrequencyAggregator, MidpointCalculator, DistanceFilter
from src.statistics.inter_statistics import InterStatisticsCalculator
from src.statistics.spline_fit import SplineDataPreparer, SplineFitter, SplineAnalysis
from src.statistics.pvals import IntraPValueCalculator, InterPValueCalculator
from src.statistics.fdr import FDRCalculator
from src.setup.config_loader import Config

class StatController:

    def __init__(self, config: Config, filtering_output) -> None:
        self.config = config
        self.metadata = filtering_output.metadata
        self.data = filtering_output.data

    def run(self):
        if self.metadata.interaction_type == "intra":
            return self._process_intra_stats()
        elif self.metadata.interaction_type == "inter":
            return self._process_inter_stats()
        else:
            raise ValueError("Invalid interaction type")

    def _process_intra_stats(self):
        filtered_data = self._prepare_intra_data(intra_data=self.data)
        spline_input = self._create_spline_input(filtered_data)
        spline = self._fit_spline(spline_input)
        pvals = self._calculate_intra_pvals(spline, spline_input)
        fdr = self._calculate_fdr_intra(pvals)
        return StatisticalOutput(metadata=self.metadata, data=fdr)

    def _process_inter_stats(self):
        inter_stats = self._prepare_inter_data(self.data)
        pvals = self._process_pvals_inter(inter_stats)
        fdr = self._process_fdr_inter(pvals)
        return StatisticalOutput(metadata=self.metadata, data=fdr)

    def _prepare_intra_data(self, intra_data):
        midpoints = MidpointCalculator(self.config).calculate_midpoints(intra_data)
        distances = DistanceCalculator(self.config).calculate_distances(midpoints)
        binned_data = EqualOccupancyBinner(self.config, intra_data).bin_data(distances)
        aggregated_data = FrequencyAggregator(self.config).aggregate_frequencies(binned_data)
        filtered_data = DistanceFilter(self.config, aggregated_data, self.metadata.resolution).filter_distances()
        return filtered_data

    def _prepare_inter_data(self, inter_data):
        inter_stats = InterStatisticsCalculator(inter_data, self.config).compute_inter_stats()
        return inter_stats

    @staticmethod
    def _create_spline_input(filtered_data):
        spline_input = SplineDataPreparer.prepare_data_for_spline(filtered_data)
        return spline_input

    def _fit_spline(self, spline_input):
        spline_fitter = SplineFitter(self.config, spline_input[0], spline_input[1], spline_input[2])
        spline = spline_fitter.fit_spline()
        return spline

    def _spline_stats(self, spline, spline_input):
        spline_analysis = SplineAnalysis(self.config, spline, spline_input)
        spline_stats = spline_analysis.calculate_mse(), spline_analysis.calculate_residuals(), spline_analysis.calculate_r_squared()
        return spline_stats

    def _calculate_intra_pvals(self, spline, spline_input):
        pval_calculator = IntraPValueCalculator(spline_input, spline, self.metadata, self.metadata.max_possible_interaction_count_intra, self.config)
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
