# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import BedpeOutput, GroupedFiles
import dask.dataframe as dd


class BedpeCreator:

    def __init__(self, config: Config, grouped_files: GroupedFiles):
        self.config = config
        self.grouped_files = grouped_files

    def run(self):
        bedpe_ddf = self._make_bedpe_dask()

        if self.config.pipeline_settings.use_hicpro_bias:
            bias_series = self._read_and_process_bias()
            bedpe_ddf = self._merge_bedpe_and_bias(bedpe_ddf, bias_series)

        interaction_calculator = InteractionCalculator(self.config, bedpe_ddf)
        total_intra, total_inter, interaction_count_per_chromosome = interaction_calculator.calculate_total_interactions()

        # Use Dask's compute method to evaluate all delayed values
        total_intra, total_inter, interaction_count_per_chromosome = dd.compute(total_intra, total_inter, interaction_count_per_chromosome)

        # Update the metadata with the interaction counts
        self.grouped_files.metadata.total_interaction_count_intra = total_intra
        self.grouped_files.metadata.total_interaction_count_inter = total_inter
        self.grouped_files.metadata.interaction_count_per_chromosome = interaction_count_per_chromosome

        return BedpeOutput(metadata=self.grouped_files.metadata, data=bedpe_ddf)

    def _make_bedpe_dask(self) -> dd.DataFrame:

        bed_dtypes = {'chr': str, 'start': 'int32', 'end': 'int32', 'idx': 'int32'}
        matrix_dtypes = {'id1': 'int32', 'id2': 'int32', 'interaction_count': 'int32'}

        bed_ddf = dd.read_csv(self.grouped_files.bed_file, sep="\t", header=None, names=["chr", "start", "end", "idx"], dtype=bed_dtypes)
        matrix_ddf = dd.read_csv(self.grouped_files.matrix_file, sep="\t", header=None, names=["id1", "id2", "interaction_count"], dtype=matrix_dtypes)

        merged_ddf = matrix_ddf.merge(bed_ddf, left_on="id1", right_on="idx", how="left") \
            .merge(bed_ddf, left_on="id2", right_on="idx", how="left", suffixes=("_1", "_2"))

        bedpe_ddf = merged_ddf[["chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2", "interaction_count"]]

        return bedpe_ddf

    def _read_and_process_bias(self):
        # Read bias file into Dask DataFrame
        bias_series = dd.read_csv(self.grouped_files.metadata.bias_file_path, header=None, squeeze=True)

        # Replace "nan" string with -1
        bias_series = bias_series.mask(bias_series == "nan", -1).astype(float)
        lower_bound = self.config.pipeline_settings.bias_lower_bound
        upper_bound = self.config.pipeline_settings.bias_upper_bound

        # Filter out-of-bounds values
        bias_series = bias_series.where(bias_series.between(lower_bound, upper_bound), -1)
        return bias_series

    @staticmethod
    def _merge_bedpe_and_bias(merged_ddf, bias_series):
        # Check if extending the bias series is necessary
        if len(bias_series) < (merged_ddf["id1"].max() + 1):

            # TODO: Log this later
            # Extend bias_series to match the length of bedpe_ddf (BED file length)
            extended_bias_series = bias_series.reindex(range(merged_ddf['id1'].max() + 1), fill_value=-1)
        else:
            # Use original bias series if it covers all IDs
            extended_bias_series = bias_series

        # Merge extended bias values with bedpe_ddf
        merged_ddf["bias_1"] = extended_bias_series.loc[merged_ddf["id1"]].values
        merged_ddf["bias_2"] = extended_bias_series.loc[merged_ddf["id2"]].values

        return merged_ddf

class InteractionCalculator:
    def __init__(self, config, data):
        self.config = config
        self.data = data  # This is the dask dataframe with interactions

    def calculate_total_interactions(self):
        chrom_counts = self.data['chr'].value_counts().compute()
        total_bins = chrom_counts.sum()

        intra_counts = (chrom_counts * (chrom_counts + 1)) // 2
        inter_counts = chrom_counts * (total_bins - chrom_counts)

        total_possible_intra = intra_counts.sum()
        total_possible_inter = inter_counts.sum() // 2

        return total_possible_intra, total_possible_inter, intra_counts.to_dict()
