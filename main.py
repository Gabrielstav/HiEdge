# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import argparse
from src.setup.config_loader import Config
from src.setup.pipeline_input import HicProInputFilePreparer
from src.setup.setup_tool import RunDirectorySetup
from src.setup.config_loader import InstantiateConfig
from src.data_preparation.preparation_controller import DataPreparationController
from src.filtering.filtering_controller import FilteringController
from src.statistics.stat_controller import StatController
from src.output.output_formatter import OutputConfiguratorRunner
from pathlib import Path
from dask import compute
from dask.delayed import delayed
from typing import List
from time import time


class Pipeline:

    def __init__(self, config_path: Path):
        self.config_path = config_path
        self.config = self._load_config()

    def run(self):

        # start time of the pipeline
        start_time = time()

        # Set up the run directory and instantiate config
        setup_tool = RunDirectorySetup(config=self.config, config_path=self.config_path)
        setup_tool.prepare_run_environment()

        # Prepare input files and create metadata
        input_file_preparer = HicProInputFilePreparer(self.config)
        grouped_files_list = input_file_preparer.prepare_input_files()
        print(f"Preparing input files: {grouped_files_list}")

        # Create interaction data
        prepared_interactions = self._execute_in_parallel(grouped_files_list, DataPreparationController)
        print(f"Making interaction datasets: {prepared_interactions}")

        # Filter on the interactions
        filtered_interactions = self._execute_in_parallel(prepared_interactions, FilteringController)
        print(f"Filtering on interaction datasets: {filtered_interactions}")

        # Perform statistical analysis
        statistical_output = self._execute_in_parallel(filtered_interactions, StatController)
        print(f"Doing statistical modeling: {statistical_output}")

        # Configure and write output
        self._execute_in_parallel(statistical_output, OutputConfiguratorRunner)
        print(f"Writing output!")

        # print time taken to run the pipeline
        print(f"Time taken to run the pipeline: {time() - start_time} seconds.")

    def _execute_in_parallel(self, inputs, controller_class) -> List:
        # Ensure inputs is a list for uniform processing
        if not isinstance(inputs, list):
            inputs = [inputs]

        # Create delayed tasks for running each controller class instance
        delayed_runs = [delayed(controller_class(self.config, input_obj).run()) for input_obj in inputs]

        print(f"Starting parallel execution of {controller_class.__name__} with {len(inputs)} input(s).")

        if not delayed_runs:
            print("No tasks to execute.")
            return []

        # Execute tasks in parallel and collect the results
        results = compute(*delayed_runs)

        # Results are returned as a tuple, convert to list if multiple results, else return single result
        return list(results) if len(inputs) > 1 else results[0]

    def _load_config(self) -> Config:
        config_mapper = InstantiateConfig(self.config_path)
        return config_mapper.load_configuration_from_file()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the pipeline with the specified configuration.")
    parser.add_argument("--config", "-con", "-c", type=Path, required=True, help="Path to the configuration YAML file.")
    parser.add_argument("--run_name", "-r", type=str, required=False, help="Name of the run.")
    args = parser.parse_args()

    pipeline = Pipeline(args.config)
    pipeline.run()
