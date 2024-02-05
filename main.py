# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from src.setup.config_loader import Config
from src.setup.data_structures import Metadata, GroupedFiles, InteractionOutput, FilteringOutput, SplineInput, StatisticalOutput
from src.data_preparation.preparation_controller import DataPreparationController
from src.filtering.filtering_controller import FilteringController
from src.statistics.stat_controller import StatController
from src.output.output_formatter import OutputConfigurator


class Pipeline:
    def __init__(self, config: Config):
        self.config = config

    def run(self):
        # Load input data and metadata
        initial_data, initial_metadata = self._load_initial_data()

        # Run data preparation
        data_prep_controller = DataPreparationController(self.config, initial_data, initial_metadata)
        data_prep_output = data_prep_controller.run()

        # Run filtering
        filtering_controller = FilteringController(self.config, data_prep_output.data, data_prep_output.metadata)
        filtering_output = filtering_controller.run_filters()

        # Run statistical analysis
        stat_controller = StatController(self.config, filtering_output.data, filtering_output.metadata)
        stat_output = stat_controller.run_stats()

        # Configure and write output
        output_configurator = OutputConfigurator(self.config)
        output_configurator.configure_output(stat_output)

    def _load_initial_data(self):
        # Replace with your actual data and metadata loading logic
        initial_data = loadData()  # Function to load your data
        initial_metadata = loadMetadata()  # Function to load your metadata
        return initial_data, initial_metadata

    def _run_setup(self):
        default_config_path = Path(__file__).parent / "default_config.yaml"
        run_directory_setup = RunDirectorySetup(default_config_path=default_config_path)

        # Replace with desired path and run name or fetch them from the config
        run_directory = run_directory_setup.create_run_directory(path=self.config.run_path, run_name=self.config.run_name)

        if run_directory:
            print(f"Run directory setup complete: {run_directory}")
        else:
            print("Run directory setup failed.")
            return

if __name__ == "__main__":
    config = loadConfig()  # Replace with your config loading logic
    pipeline = Pipeline(config)
    pipeline.run()




