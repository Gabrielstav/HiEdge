# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import argparse
from src.setup.config_loader import Config
from src.setup.pipeline_input import HicProInputFilePreparer
from src.setup.data_structures import Metadata, GroupedFiles, InteractionOutput, FilteringOutput, SplineInput, StatisticalOutput
from src.setup.setup_tool import RunDirectorySetup
from src.setup.config_loader import ConfigMapper
from src.data_preparation.preparation_controller import DataPreparationController
from src.filtering.filtering_controller import FilteringController
from src.statistics.stat_controller import StatController
from src.output.output_formatter import OutputConfigurator
from pathlib import Path


class Pipeline:

    def __init__(self, config_path: Path):
        self.config_path = config_path
        self.config = self._load_config()

    def run(self):
        # Set up the run directory
        setup_tool = RunDirectorySetup(config=self.config)
        setup_tool.prepare_run_environment()

        # Prepare Hi-C Pro input files
        input_file_preparer = HicProInputFilePreparer(self.config)
        grouped_files_list = input_file_preparer.prepare_input_files()

        # Data preparation
        data_prep_controller = DataPreparationController(self.config, grouped_files_list)
        data_prep_output = data_prep_controller.run()

        # Filtering
        filtering_controller = FilteringController(self.config, data_prep_output.data, data_prep_output.metadata)
        filtering_output = filtering_controller.run_filters()

        # Statistical analysis
        stat_controller = StatController(self.config, filtering_output.data, filtering_output.metadata)
        stat_output = stat_controller.run_stats()

        # Configure and write output
        output_configurator = OutputConfigurator(self.config)
        output_configurator.configure_output(stat_output)

        # Print final message with runtime and output directory
        print(f"Pipeline finished. Output files written to {self.config.paths.output_dir}")

    def _load_config(self) -> Config:
        config_mapper = ConfigMapper(self.config_path)
        return config_mapper.load_configuration_from_file()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the pipeline with the specified configuration.")
    parser.add_argument("--config", type=Path, required=True, help="Path to the configuration YAML file.")
    args = parser.parse_args()

    pipeline = Pipeline(args.config)
    pipeline.run()
