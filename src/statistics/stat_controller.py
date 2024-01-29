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

    def __init__(self, config: Config, data: dd.DataFrame, metadata: Metadata):
        self.config = config
        self.data = data
        self.metadata = metadata

    def run_statistical_tests(self):
        # Depending on interaction type, run different statistical tests
        if self.metadata.interaction_type == "intra":
            return self._run_intra_tests()
        elif self.metadata.interaction_type == "inter":
            return self._run_inter_tests()
        else:
            raise ValueError("Unknown interaction type")

    def _run_intra_tests(self):
        # Step 1: Compute intra statistics and prepare data for spline fitting
        intra_statistics = IntraStatistics(self.config, self.data)
        prepared_data = intra_statistics.prepare_data()

        # Step 2: Fit spline to the prepared data
        spline_fitter = SplineFitter(prepared_data, self.config)
        spline = spline_fitter.fit()

        # Step 3: Calculate p-values using the fitted spline
        p_value_calculator = IntraPValueCalculator(prepared_data, spline, self.metadata.max_possible_interaction_count_intra, self.config)
        p_value_results = p_value_calculator.run()

        # Step 4: Calculate FDR
        fdr_calculator = FDRCalculator(self.config, p_value_results)
        fdr_results = fdr_calculator.run()

        # Step 5: Combine results and return
        return StatisticsOutput(metadata=self.metadata, data=fdr_results)

    def _run_inter_tests(self):
        # Step 1: Compute inter statistics
        inter_statistics = InterStatistics(self.config, self.data)
        statistics_results = inter_statistics.compute_statistics()

        # Step 2: Calculate p-values for interchromosomal data
        p_value_calculator = InterPValueCalculator(self.config, statistics_results, self.metadata.max_possible__interaction_count_inter)
        p_value_results = p_value_calculator.run()

        # Step 3: Calculate FDR
        fdr_calculator = FDRCalculator(self.config, p_value_results)
        fdr_results = fdr_calculator.run()

        # Step 4: Combine results and return
        return StatisticsOutput(metadata=self.metadata, data=fdr_results)

# Data class for statistical test outputs


