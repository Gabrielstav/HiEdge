# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import shutil
import datetime
from pathlib import Path
from dataclasses import asdict
import yaml
from src.setup.config_loader import Config


class RunDirectorySetup:
    subdirs = ["output", "logs", "tmp", "config"]

    def __init__(self, config: Config):
        self.config = config
        # Ensure run_name is incorporated into the final run directory path
        self.run_name = self.config.run_name or f"run_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.run_path = self.config.paths.run_dir / self.run_name  # Adjusted to include self.run_name in the path

    def prepare_run_environment(self):
        if self.run_path.exists():
            print(f"Overwriting files in directory: {self.run_path}")
            shutil.rmtree(self.run_path)

        self.run_path.mkdir(parents=True)
        for subdir in self.subdirs:
            (self.run_path / subdir).mkdir()

        if not all((self.config.paths.run_dir / subdir).exists() for subdir in self.subdirs):
            print("Error in setting up run directory. Not all subdirectories were created.")

        print(f"Run directory created at {self.run_path}")
        self.copy_config_file()

    def copy_config_file(self):
        config_filename = f"config_{self.run_name}.yaml"
        config_target_path = self.run_path / "config" / config_filename
        with config_target_path.open("w") as file:
            # Serialize the config dataclass to a dictionary before dumping
            yaml.safe_dump(asdict(self.config), file, default_flow_style=False)
        print(f"Configuration file used in current run copied to {config_target_path}")
