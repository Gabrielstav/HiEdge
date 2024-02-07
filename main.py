# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import argparse
from src.setup.config_loader import Config
from src.setup.pipeline_input import HicProInputFilePreparer
from src.setup.setup_tool import RunDirectorySetup
from src.setup.config_loader import ConfigMapper
from src.data_preparation.preparation_controller import DataPreparationController
from src.filtering.filtering_controller import FilteringController
from src.statistics.stat_controller import StatController
from src.output.output_formatter import OutputConfigurator
from pathlib import Path
from dask import delayed, compute
from typing import List


class Pipeline:

    def __init__(self, config_path: Path):
        self.config_path = config_path
        self.config = self._load_config()

    def run(self):

        # Set up the run directory and instantiate config
        setup_tool = RunDirectorySetup(config=self.config)
        setup_tool.prepare_run_environment()

        # Prepare input files and create metadata
        input_file_preparer = HicProInputFilePreparer(self.config)
        grouped_files_list = input_file_preparer.prepare_input_files()

        # Create interaction data
        prepared_interactions = self._execute_in_parallel(grouped_files_list, DataPreparationController)

        # Filter on the interactions
        filtered_interactions = self._execute_in_parallel(prepared_interactions, FilteringController)

        # Perform statistical analysis
        statistical_output = self._execute_in_parallel(filtered_interactions, StatController)

        # Configure and write output
        self._execute_in_parallel(statistical_output, OutputConfigurator)

    def _execute_in_parallel(self, inputs, controller_class) -> List:
        # Create delayed objects for each contrller class
        delayed_tasks = [delayed(controller_class)(self.config, input_obj) for input_obj in inputs]
        delayed_runs = [task.run() for task in delayed_tasks]

        # Execute in parallel and return the results
        return compute(*delayed_runs)[0]  # returns the first element of compute tuple, list of results

    def _load_config(self) -> Config:
        config_mapper = ConfigMapper(self.config_path)
        return config_mapper.load_configuration_from_file()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the pipeline with the specified configuration.")
    parser.add_argument("--config", type=Path, required=True, help="Path to the configuration YAML file.")
    args = parser.parse_args()

    pipeline = Pipeline(args.config)
    pipeline.run()
