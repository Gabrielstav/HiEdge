# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
import dask.dataframe as dd

class FilterChromosomes:

    def __init__(self, config: Config):
        self.config = config
        self.chromosome_action = self._determine_chromosome_action()
        self.chromosomes = self._get_chromosomes()

    def _determine_chromosome_action(self):
        if self.config.pipeline_settings.select_chromosomes:
            return "select"
        elif self.config.pipeline_settings.remove_chromosomes:
            return "remove"
        else:
            return None

    def _get_chromosomes(self):
        return self.config.pipeline_settings.select_chromosomes if self.chromosome_action == "select" \
            else self.config.pipeline_settings.remove_chromosomes

    def filter_data(self, data: dd.DataFrame) -> dd.DataFrame:
        if self.chromosome_action == "select":
            return self._filter_select_chromosomes(data)
        elif self.chromosome_action == "remove":
            return self._filter_remove_chromosomes(data)
        return data

    def _filter_select_chromosomes(self, data: dd.DataFrame) -> dd.DataFrame:
        is_intra = data["interaction_type"] == "intra"
        is_inter = data["interaction_type"] == "inter"
        in_chromosomes = data["chr_1"].isin(self.chromosomes) | data["chr_2"].isin(self.chromosomes)

        return data[(is_intra & (data["chr_1"].isin(self.chromosomes))) |
                    (is_inter & in_chromosomes)]

    def _filter_remove_chromosomes(self, data: dd.DataFrame) -> dd.DataFrame:
        is_intra = data["interaction_type"] == "intra"
        is_inter = data["interaction_type"] == "inter"
        not_in_chromosomes = ~data["chr_1"].isin(self.chromosomes) & ~data["chr_2"].isin(self.chromosomes)

        return data[(is_intra & (~data["chr_1"].isin(self.chromosomes))) |
                    (is_inter & not_in_chromosomes)]
