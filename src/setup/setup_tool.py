# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import shutil
import datetime
from src.setup.config_loader import Config
from pathlib import Path


class RunDirectorySetup:
    subdirs = ["output", "logs", "tmp", "config", "plots"]

    def __init__(self, config: Config, config_path: Path):
        self.config = config
        self.config_path = config_path
        self.run_name = self.config.run_name or f"run_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.run_path = self.config.paths.run_dir / self.run_name

    def prepare_run_environment(self):
        if self.run_path.exists():
            print(f"Overwriting files in directory: {self.run_path}")
            shutil.rmtree(self.run_path)

        self.run_path.mkdir(parents=True)
        for subdir in self.subdirs:
            (self.run_path / subdir).mkdir()

        if not all((self.run_path / subdir).exists() for subdir in self.subdirs):
            print(f"Error in setting up run directory. Not all subdirectories were created: {self.subdirs}")

        print(f"Run directory created at {self.run_path}")
        self.copy_config_file()

        # update the config instance with the correct run path
        self.config.paths.run_dir = self.run_path

    def copy_config_file(self):
        config_filename = f"config_{self.run_name}.yaml"
        config_target_path = self.run_path / "config" / config_filename
        shutil.copy(self.config_path, config_target_path)
        print(f"Configuration file used in current run copied to {config_target_path}")
