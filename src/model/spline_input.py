# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import ConfigMapper as Config
from src.setup.data_structures import FilteringOutput, BlacklistOutput, CytobandOutput, StatInput

"""
DataProcessor Class: Handles loading and organizing Hi-C data.
ProbabilityCalculator Class: Responsible for computing interaction probabilities, calculate the mean interaction frequency per bin and prepare data for spline fitting.
SplineFitter Class: Encapsulates the spline fitting process. This class will take the output from ProbabilityCalculator and fit a spline to model the expected interaction frequencies as a function of genomic distance.
StatisticalModel Class: Performs statistical testing. It will use the spline model to calculate p-values for each interaction, considering the total sum of contact counts and observed count for a given interaction.
ResultAggregator Class: Gathers and compiles results, particularly for parallel computation scenarios.

"""

class ResolveBedpeInput:
    """
    Resolves what BEDPE data class to pass to processing for statistical modeling based on config file (FilteringOutput, BlacklistOutput, CytobandOutput)
    using the config instance to decide.
    :return: Dataclass containing data and metadata
    """
    def __init__(self, config: Config, data_output):
        self.config = config
        self.data_output = data_output

    def resolve_input(self):
        if isinstance(self.data_output, CytobandOutput) and self.config.config_data.pipeline_settings.filter_cytobands:
            return self.data_output
        elif isinstance(self.data_output, BlacklistOutput) and self.config.config_data.pipeline_settings.filter_blacklist and not isinstance(self.data_output, CytobandOutput):
            return self.data_output
        elif isinstance(self.data_output, FilteringOutput):
            return self.data_output
        else:
            raise ValueError("Invalid data class type provided")

    # Need to figure out what data to process (depending on what filtering classes are used) and return the appropriate dataclass instance contaning the data.
    # We either use the data from the FilteringOutput class or the BlacklistOutput class or the CytobandOutput class.

class CalculateMidpoints:

    def __init__(self, config: Config, data_input):
        self.config = config
        self.data_ddf = data_input.bedpe_ddf
        self.metadata = data_input.metadata

    def process_data(self):
        # Calculate midpoints for each bin using vectorized operations if metadata == intra
        self.data_ddf['midpoint_1'] = (self.data_ddf['start_1'] + self.data_ddf['end_1']) // 2
        self.data_ddf['midpoint_2'] = (self.data_ddf['start_2'] + self.data_ddf['end_2']) // 2

        # Create and return StatInput dataclass instance with the updated DataFrame
        stat_input = StatInput(metadata=self.metadata, bedpe_ddf=self.data_ddf)
        return stat_input

    pass

class CalculateGenomicDistance:

        def __init__(self, config: Config, data_input):
            self.config = config
            self.data_ddf = data_input.bedpe_ddf
            self.metadata = data_input.metadata

        pass

class MeanInteractionFrequency:

    def __init__(self, config: Config, data_input):
            self.config = config
            self.data_ddf = data_input.bedpe_ddf
            self.metadata = data_input.metadata

    pass

class ApplyBiasToCounts:

    def __init__(self, config: Config, data_input):
            self.config = config
            self.data_ddf = data_input.bedpe_ddf
            self.metadata = data_input.metadata

    pass

