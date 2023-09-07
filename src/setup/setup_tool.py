# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import shutil
import datetime
import argparse
from pathlib import Path
from typing import Optional
import yaml

class RunDirectorySetup:

    SUBDIRS = ["output", "logs", "tmp", "config"]

    def __init__(self, default_config_path: Path):
        self.default_config_path = default_config_path

    @staticmethod
    def adjust_config_paths(config_path: Path, run_directory: Path):
        """
        Adjust the paths in the configuration file according to the run directory.
        """
        try:
            with config_path.open('r') as file:
                config = yaml.safe_load(file)

            # Adjust paths
            for key, value in config['paths'].items():
                if value.startswith("./"):  # It's a relative path
                    config['paths'][key] = str(run_directory / Path(value).name)

            with config_path.open('w') as file:
                yaml.safe_dump(config, file, default_flow_style=False)

        except Exception as e:
            print(f"Error processing config file: {e}")
            return

    def create_run_directory(self, path: Optional[Path] = None, run_name: Optional[str] = None) -> Optional[Path]:
        """
        Creates a run directory with a given name at a specified path.
        If no name or path is provided, it uses defaults.
        """
        # If no path is provided, use current working directory
        if not path:
            path = Path.cwd()

        # If no run_name is provided, use timestamp as name
        if not run_name:
            timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            run_name = f"run_{timestamp}"

        run_path = path / run_name
        config_subdir = run_path / "config"

        # Check if the dir for the run exists
        if run_path.exists():
            user_input = input(f"Directory {run_path} already exists! Do you wish to overwrite? (y/n): ")
            if user_input.lower() == 'y':
                shutil.rmtree(run_path)
            else:
                print("Setup halted.")
                return

        # Create the run dir and standard sub-dirs
        run_path.mkdir(parents=True)
        for subdir in ["output", "logs", "tmp", "config"]:
            (run_path / subdir).mkdir()

        # Copy the default configuration file to the run directory's config sub-dir
        config_filename = f"config_{run_name}.yaml"
        shutil.copy(self.default_config_path, config_subdir / config_filename)

        self.adjust_config_paths(config_subdir / config_filename, run_path)

        print(f'Run configuration "{run_name}" created in "{run_path}"')
        print(f"You can modify the config file at: {config_subdir / config_filename}")

        return run_path


def main():
    parser = argparse.ArgumentParser(description="Set up a new run directory.")

    # Args for the CLI
    parser.add_argument('-p', '--path', type=Path, help='Path where the run directory should be created. Default is current directory.')
    parser.add_argument('-n', '--name', type=str, help='Name for the run directory. Default is "run_{timestamp}".')

    args = parser.parse_args()

    # Ensure that the default config file is present.
    default_config_path = Path(__file__).parent / 'default_config.yaml'
    if not default_config_path.exists():
        print(f"Error: Default configuration file not found at {default_config_path}. Please ensure it's in the correct location.")
        return

    # Creating the run directory using the SetUpTool
    tool = RunDirectorySetup(default_config_path=default_config_path)
    tool.create_run_directory(path=args.path, run_name=args.name)


if __name__ == "__main__":
    main()

