# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import SplineInput, FilteringOutput, Metadata
from src.statistics.intra_statistics import DataPreparationController
from src.statistics.inter_statistics import DataPreparationController
from src.statistics.fdr import FDRCalculator
from src.setup.data_structures import FilteringOutput
from src.setup.config_loader import Config

class StatisticalTestController:

    def __init__(self, config: Config, data):
        self.config = config
        self.data = data
        self.p_value_calculator = IntraPValueCalculator(config, data)
        self.fdr_calculator = None  # To be instantiated later with p-value results

    def run_statistical_tests(self):

        self.data = self.p_value_calculator.run()

        self.fdr_calculator = FDRCalculator(self.config, self.data)
        fdr_results = self.fdr_calculator.run()

        # Combine p-value and FDR results
        combined_results = self.combine_results(self.data, fdr_results)

        return combined_results

    @staticmethod
    def combine_results(data, fdr_results):
        # Assuming both dataframes are aligned and have the same index
        combined_data = data.copy()
        combined_data["adjusted_p_value"] = fdr_results["adjusted_p"]
        return combined_data

class StatisticsController:

    def __init__(self, config, data, interaction_type):
        self.config = config
        self.data = data
        self.interaction_type = interaction_type

    def run(self):
        if self.interaction_type == "intra":
            return self._process_intra_data()
        elif self.interaction_type == "inter":
            return self._process_inter_data()

    def _process_intra_data(self):
        intra_stats = IntraStatisticsCalculator(self.config, self.data).calculate_statistics()
        spline = SplineFitter(intra_stats).fit_spline()
        p_values = IntraPValueCalculator(intra_stats, spline).calculate_p_values()
        fdr_results = IntraFDRCalculator(self.config, p_values).calculate_fdr()
        return fdr_results

    def _process_inter_data(self):
        inter_stats = InterStatisticsCalculator(self.config, self.data).calculate_statistics()
        p_values = InterPValueCalculator(inter_stats).calculate_p_values()
        fdr_results = InterFDRCalculator(self.config, p_values).calculate_fdr()
        return fdr_results
