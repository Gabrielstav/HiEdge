# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.data_structures import SplineInput, FilteringOutput, Metadata, StatisticalOutput
from src.statistics.intra_statistics import DistanceCalculator, EqualOccupancyBinner, FrequencyAggregator, MidpointCalculator, DistanceFilter
from src.statistics.inter_statistics import InterStatisticsCalculator
from src.statistics.spline_fit import SplineDataPreparer, SplineFitter, SplineAnalysis
from src.statistics.pvals import IntraPValueCalculator, InterPValueCalculator
from src.statistics.fdr import FDRCalculator
from src.setup.config_loader import Config
import dask.dataframe as dd

class StatController:

    def __init__(self, config: Config, data: dd.DataFrame, metadata: Metadata):
        self.config = config
        self.data = data
        self.metadata = metadata

    def run_stats(self):
        # Call relevant methods in sequence to run all stats

        return StatisticalOutput(metadata=self.metadata, data=self.data)

    def _process_intra_data(self, intra_data):
        pass

    def _process_inter_data(self, inter_data):
        pass

    def _process_spline(self, spline_input):
        pass

    def _process_pvals(self, pval_input):
        pass

    def _process_fdr(self, fdr_input):
        pass
