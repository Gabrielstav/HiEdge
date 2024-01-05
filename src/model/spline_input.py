# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import BlacklistOutput, CytobandOutput, SplineInput
import dask.dataframe as dd

class ResolveBedpeInput:

    def __init__(self, config: Config, data_output):
        self.config = config
        self.data_output = data_output

    def resolve_input(self):
        # Mapping configuration options to data classes
        data_class_mapping = {
            "cytoband": (CytobandOutput, self.config.pipeline_settings.filter_cytobands),
            "blacklist": (BlacklistOutput, self.config.pipeline_settings.filter_blacklist)
        }

        for data_class, (output_type, use_filter) in data_class_mapping.items():
            if isinstance(self.data_output, output_type) and use_filter:
                return output_type

        # Return the data field of the resolved data class
        return self.data_output.data

class DataPreparationController:

    def __init__(self, config: Config, input_data):
        self.config = config
        self.input_data = input_data

    def run(self):
        # Resolve the input data class based on configuration
        input_resolver = ResolveBedpeInput(self.config, self.input_data)
        resolved_data = input_resolver.resolve_input()

        # Prepare data for spline fitting
        prepared_data = self._prepare_data_for_spline(resolved_data)

        # Trigger computation and instantiate dataclass with results
        if isinstance(prepared_data, dd.DataFrame):
            computed_data = prepared_data.compute()
        else:
            computed_data = prepared_data

        # Instantiate and return the dataclass with computed data and original metadata
        return SplineInput(metadata=self.input_data.metadata, data=computed_data)

    def _prepare_data_for_spline(self, data) -> dd.DataFrame:
        midpoint_calculator = MidpointCalculator(self.config)
        data_with_midpoints = midpoint_calculator.calculate_midpoints(data)

        distance_calculator = DistanceCalculator(self.config)
        data_with_distances = distance_calculator.calculate_distances(data_with_midpoints)

        binner = EqualOccupancyBinner(self.config, data_with_distances)
        binned_data = binner.bin_data(data_with_distances)

        aggregator = FrequencyAggregator(self.config)
        aggregated_data = aggregator.aggregate_frequencies(binned_data)

        return aggregated_data


class MidpointCalculator:

    def __init__(self, config: Config):
        self.config = config

    @staticmethod
    def calculate_midpoints(data: dd.DataFrame) -> dd.DataFrame:
        data["midpoint_1"] = ((data["start_1"] + data["end_1"]) / 2).astype("int32")
        data["midpoint_2"] = ((data["start_2"] + data["end_2"]) / 2).astype("int32")
        return data.persist()

class DistanceCalculator:

    def __init__(self, config: Config):
        self.config = config

    @staticmethod
    def calculate_distances(data: dd.DataFrame) -> dd.DataFrame:
        data["genomic_distance"] = abs(data["midpoint_1"] - data["midpoint_2"]).astype("int32")
        return data.persist()

class EqualOccupancyBinner:

    def __init__(self, config, input_data):
        self.config = config
        self.input_data = input_data

    @staticmethod
    def sort_data_by_distance(input_data) -> dd.DataFrame:
        return input_data.sort_values("genomic_distance")

    @staticmethod
    def calculate_total_contacts(data) -> dd.DataFrame:
        return data["interaction_count"].sum().compute()

    def assign_to_bins(self, sorted_data, total_contacts):
        target_per_bin = total_contacts / self.config.number_of_bins
        sorted_data["cumulative_count"] = sorted_data["interaction_count"].cumsum()
        sorted_data["bin_label"] = (sorted_data["cumulative_count"] / target_per_bin).astype("int")

        # Handle the tiebreak condition
        sorted_data = self.apply_tiebreak_condition

        return sorted_data

    def apply_tiebreak_condition(self):
        # Identify rows where the tiebreak condition applies
        self.input_data['needs_tiebreak'] = (
                (self.input_data["genomic_distance"] == self.input_data["genomic_distance"].shift(1)) &
                (self.input_data["bin_label"] != self.input_data["bin_label"].shift(1))
        )

        # Vectorized approach to adjust bin labels
        # The logic here is to 'push forward' the bin label when a tiebreak is needed
        self.input_data["adjusted_bin_label"] = self.input_data["bin_label"].where(~self.input_data["needs_tiebreak"], self.input_data["bin_label"].shift(1))

        # Clean up and finalize
        self.input_data = self.input_data.drop(["needs_tiebreak"], axis=1)
        self.input_data = self.input_data.rename(columns={"adjusted_bin_label": "bin_label"})

        return self

    @staticmethod
    def calculate_bin_statistics(binned_data):
        # Calculate average genomic distance and contact probability for each bin
        bin_stats = binned_data.groupby("bin_label").agg({
            "genomic_distance": "mean",
            "interaction_count": ["mean", "sum"]
        }).compute()
        bin_stats.columns = ["average_distance", "average_contact_probability", "total_contacts"]
        return bin_stats

    def bin_data(self, data):
        sorted_data = self.sort_data_by_distance(data)
        total_contacts = self.calculate_total_contacts(data)
        binned_data = self.assign_to_bins(sorted_data, total_contacts)
        bin_stats = self.calculate_bin_statistics(binned_data)
        return bin_stats

class FrequencyAggregator:

    def __init__(self, config: Config):
        self.config = config

    # Aggregate frequencies of genomic distances
    @staticmethod
    def aggregate_frequencies(data: dd.DataFrame) -> dd.DataFrame:
        aggregated_data = data.groupby("genomic_distance").agg({
            "interaction_count": ["sum", "mean", "std"]
        }).reset_index()
        aggregated_data.columns = ["genomic_distance", "interaction_sum", "interaction_mean", "interaction_std"]
        return aggregated_data.persist()



