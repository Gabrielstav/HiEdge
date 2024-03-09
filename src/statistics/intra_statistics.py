# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
import dask.dataframe as dd
import numpy as np
import pandas as pd


# TODO:
#   So the metabins is divided into equal occupancy bins across the genome, but the target should insteaed be that each chromosome has the same number of metabins?
#   Because in low resolution, smaller chromosomes will have all interactions in one metabin
#   Another approach is find metabin with lowest interaction count, then adjust all other metabins to have the similar interaction count
#   This will ensure that the metabins are more evenly distributed across the genome
#   But chromosomes are treated as separate entities in intra-datasets, so this might not be the best approach, as smaller chromosomes will then dictate the binning of larger ones
#

class MidpointCalculator:

    def __init__(self, config: Config):
        self.config = config

    @staticmethod
    def calculate_midpoints(data: dd.DataFrame) -> dd.DataFrame:
        data["midpoint_1"] = ((data["start_1"] + data["end_1"]) / 2).astype("int32")
        data["midpoint_2"] = ((data["start_2"] + data["end_2"]) / 2).astype("int32")
        return data


class DistanceCalculator:

    def __init__(self, config: Config):
        self.config = config

    @staticmethod
    def calculate_distances(data: dd.DataFrame) -> dd.DataFrame:
        # Calculate the genomic distance
        data["genomic_distance"] = abs(data["midpoint_1"] - data["midpoint_2"]).astype("int32")

        # Filter out data with genomic distance == 0 (self-interactions)
        filtered_data = data[data["genomic_distance"] != 0]

        return filtered_data


class DistanceFilter:

    def __init__(self, config: Config, data: dd.DataFrame, resolution: int):
        self.config = config
        self.data = data
        self.resolution = resolution

    def filter_distances(self) -> dd.DataFrame:
        if self.resolution in self.config.pipeline_settings.interaction_distance_filters and self.config.pipeline_settings.use_interaction_distance_filters:
            distance_filter_config = self.config.pipeline_settings.interaction_distance_filters[self.resolution]
            lower_bound = distance_filter_config.min_distance
            upper_bound = distance_filter_config.max_distance
            self.data = self.data[(self.data["genomic_distance"] >= lower_bound) & (self.data["genomic_distance"] <= upper_bound)]
        return self.data


class EqualOccupancyBinner:

    def __init__(self, config: Config, input_data: dd.DataFrame):
        self.config = config
        self.input_data = input_data

    def bin_data(self, data: dd.DataFrame, total_interactions: int, intra_counts_per_chromosome: dict) -> dd.DataFrame:
        sorted_data = self._sort_data_by_distance(data)
        total_contacts = self._calculate_total_contacts(data)
        binned_data = self._assign_to_bins(sorted_data, total_contacts, self.config.statistical_settings.metabin_occupancy)
        bin_stats = self._calculate_metabin_stats(binned_data, total_interactions, intra_counts_per_chromosome)

        print(f"Bin stats: {bin_stats}")
        print(f"Data columns: {data.columns}")
        return bin_stats

    @staticmethod
    def _sort_data_by_distance(input_data: dd.DataFrame) -> dd.DataFrame:
        return input_data.sort_values("genomic_distance")

    @staticmethod
    def _calculate_total_contacts(data: dd.DataFrame) -> int:
        return data["interaction_count"].sum()

    @staticmethod
    def _assign_to_bins(sorted_data: dd.DataFrame, total_contacts: int, num_bins: int) -> dd.DataFrame:
        target_per_bin = total_contacts / num_bins
        sorted_data["cumulative_count"] = sorted_data["interaction_count"].cumsum()
        sorted_data["metabin_label"] = (sorted_data["cumulative_count"] / target_per_bin).astype("int")

        # Handle the tiebreak condition
        sorted_data = EqualOccupancyBinner._apply_tiebreak_condition(sorted_data)

        return sorted_data

    @staticmethod
    def _apply_tiebreak_condition(data: dd.DataFrame) -> dd.DataFrame:
        # Determine if tiebreak is needed
        data["needs_tiebreak"] = (
                (data["genomic_distance"] == data["genomic_distance"].shift(1)) &
                (data["metabin_label"] != data["metabin_label"].shift(1))
        )

        # Adjust bin labels where tiebreak is needed
        condition = ~data["needs_tiebreak"]
        # Directly update metabin_label with adjusted values
        data["metabin_label"] = np.where(condition, data["metabin_label"], data["metabin_label"].shift(1).fillna(-1).astype(int))

        data = data.drop(["needs_tiebreak"], axis=1)

        if not isinstance(data, dd.DataFrame):
            data = dd.from_pandas(data, npartitions=1)

        return data

    @staticmethod
    def _calculate_unique_distances_of_total_possible_bins(data: dd.DataFrame) -> dd.DataFrame:
        # Calculate unique distances per chromosome
        unique_distances_per_chromosome = data[data["chr_1"] == data["chr_2"]].groupby("chr_1")["genomic_distance"].nunique().compute()
        print(f"Unique distances per chromosome: {unique_distances_per_chromosome}")
        return unique_distances_per_chromosome

    # TODO: Bin stats needs to calculate these values:

    # 0. range of distances in each metabin (EZ)
    #       - min and max distance for each bin (range)

    # 1. no. of possible pairs w/in this range of distances (harder but easy)
    #       - number of possible pairs for each distance in metabin
    #           - in interaction_calculator make method to find each distance for every possible pair
    #       * possPairsInRange = currBin[1]

    # 2. sumoverallContactCounts (EZ)
    #       - sum of interaction counts in each metabin

    # 3. Sumoveralldistances in this bin in distScaling vals (EZ)
    #       - sum of genomic distance in each metabin scaled by distScaling factor
    #       * distScaling=1000000.0
    #       * currBin[3]+=(float(intxnDistance/distScaling)*npairs)

    # 4. average contact count in each metabin (DONE)
    #       - average interaction count in each metabin (contact probability).
    #       * computed as total interaction counts in metabin divided by the number of possible pairs,
    #       * then normalized by some factor reflective of the overall dataset's interaction density.
    #       * sumCC = currBin[2], possPairsInRange = currBin[1]
    #       * observedIntraInRangeSum += interxn.getCount() (total interactions across all bins)
    #       * avgCC = (1.0*sumCC/possPairsInRange)/observedIntraInRangeSum

    # 5. average genomic distance in each metabin (DONE)
    #       - average genomic distance in each metabin
    #       * computed as the sum of genomic distances in each metabin divided by the number of possible pairs?
    #       * sumDistB4Scaling = currBin[3]
    #       * avgDist = distScaling*(sumDistB4Scaling/currBin[7])
    #       * currBin[5]=avgDist

    # 6. bins (EZ)
    #       - number of bins (or just the bin ID's) in each metabin?

    # TODO: Figure out if this is applicable and if so, how it works
    # 7. no. of possible pairs w/ proper dist (harder)
    #       - similar to #1, but increasing or decreasing based on possible contacts withtin the distance range?

    @staticmethod
    def _calculate_metabin_stats(data_binned: dd.DataFrame, total_interactions: int, intra_counts_per_chromosome: dict) -> dd.DataFrame:
        # Get max possible interacting bins per chromosome from metadata
        counts_df = pd.DataFrame(list(intra_counts_per_chromosome.items()), columns=["chr_1", "possible_pairs"])

        # get unique distances per chromosome
        unique_distances_per_chromosome = data_binned[data_binned["chr_1"] == data_binned["chr_2"]].groupby("chr_1")["genomic_distance"].nunique().compute()
        print(f"Unique distances per chromosome: {unique_distances_per_chromosome}")

        # TODO: Possible pairs is the same across metabins, need to adjust this based on genomic distance in metabins

        # Ensure bin_stats is a DataFrame that includes a "chr_1" column for merging
        bin_stats = data_binned.groupby("metabin_label").agg({
            "genomic_distance": "mean",
            "interaction_count": "sum",
            "chr_1": "first"  # If all entries in one metabin belong to the same chromosome
        }).reset_index().rename(columns={
            "genomic_distance": "average_distance",
            "interaction_count": "total_possible_pairs"
        })

        # Merge to associate each bin with its possible pairs count
        bin_stats = dd.merge(bin_stats, counts_df, on="chr_1", how="left")

        # Calculate average contact probability per metabin
        # TODO: This is wrong, should at least use per chromosome total interactions
        bin_stats["average_contact_probability"] = bin_stats["total_possible_pairs"] / total_interactions

        print(f"Head and tail of bin stats ddf: {bin_stats.compute()}")
        return bin_stats


# TODO: This is not called, fix this in stat_controller.py?
class FrequencyAggregator:

    def __init__(self, config: Config):
        self.config = config

    # Aggregate frequencies of genomic distances
    @staticmethod
    def aggregate_frequencies(data: dd.DataFrame) -> dd.DataFrame:
        aggregated_data = data.groupby("average_distance").agg({
            "total_possible_bins": ["sum", "mean", "std"]
        }).reset_index()
        aggregated_data.columns = ["average_distance", "interaction_sum", "interaction_mean", "interaction_std"]
        return aggregated_data
