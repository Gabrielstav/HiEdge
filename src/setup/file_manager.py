# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import os
import re
import pathlib as pl
import shutil as shutil
from typing import List, Set, Tuple, Dict, Union
from config_loader import ConfigLoader

class FileManager:

    raw_subdirectory_name = "raw"
    iced_subdirectory_name = "iced"

    @staticmethod
    def create_directory(dir_path):
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.mkdir(dir_path)

    @staticmethod
    def validate_files_exist(files):
        for file in files:
            if not os.path.isfile(file):
                print(f"File {file} does not exist.")

    @staticmethod
    def save_to_file(data, file_path):
        with open(file_path, "w") as f:
            f.writelines(data)

    @staticmethod
    def find_subdirectories_by_name(root_directory, subdirectory_name):
        return [root for root, _, _ in os.walk(root_directory) if os.path.basename(root) == subdirectory_name]

    @staticmethod
    def find_raw_subdirectories(root_directory):
        return FileManager.find_subdirectories_by_name(root_directory, FileManager.raw_subdirectory_name)

    @staticmethod
    def find_iced_subdirectories(root_directory):
        return FileManager.find_subdirectories_by_name(root_directory, FileManager.iced_subdirectory_name)

    @staticmethod
    def filter_files_on_resolution(input_files, found_resolutions_in, resolutions_provided):
        filtered_files = []
        for file in input_files:
            resolution = FileNameFormatter.extract_resolution(file)
            if resolutions_provided is None or resolution in resolutions_provided:
                filtered_files.append(file)
            found_resolutions_in.add(resolution)
        return filtered_files

    @staticmethod
    def collect_files_from_subdirectories(subdirectories, file_extension, resolutions_provided):
        found_resolutions = set()
        files = []
        for subdirectory_path in subdirectories:
            for root, _, files_in_directory in os.walk(subdirectory_path):
                relevant_files = [os.path.join(root, file) for file in files_in_directory if file.endswith(file_extension)]
                files += FileManager.filter_files_on_resolution(relevant_files, found_resolutions, resolutions_provided)
        return files, found_resolutions

    @staticmethod
    def find_files(*root_directories):
        resolutions_provided = SetDirectories.get_resolutions() # YAML config instead

        raw_subdirectories = [FileManager.find_raw_subdirectories(directory) for directory in root_directories]
        iced_subdirectories = [FileManager.find_iced_subdirectories(directory) for directory in root_directories]

        bedfiles, _ = FileManager.collect_files_from_subdirectories(raw_subdirectories, ".bed", resolutions_provided)
        matrixfiles, found_resolutions = FileManager.collect_files_from_subdirectories(raw_subdirectories, ".matrix", resolutions_provided)
        iced_matrixfiles, found_resolutions_iced = FileManager.collect_files_from_subdirectories(iced_subdirectories, ".matrix", resolutions_provided)

        FileManager.handle_missing_files(found_resolutions.union(found_resolutions_iced), resolutions_provided, len(bedfiles), len(matrixfiles), len(iced_matrixfiles))

        return bedfiles, matrixfiles, iced_matrixfiles

    @staticmethod
    def handle_missing_files(bedfile_count, matrixfile_count, iced_matrixfile_count, found_resolutions, resolutions_provided):

        # Check if no files were found for the provided resolutions
        if all(count == 0 for count in [bedfile_count, matrixfile_count, iced_matrixfile_count]) and resolutions_provided is not None:
            print(f"No files found for the provided resolutions: {resolutions_provided}. "
                  f"Resolutions found in the data: {found_resolutions}")

        # Check if no iced or raw matrix files were found corresponding to the provided command line arg
        if SetDirectories.get_normalized_data() and iced_matrixfile_count == 0:
            print(f"No ICE-normalized matrix files found, but command line argument -m set to True")
        elif not SetDirectories.get_normalized_data() and matrixfile_count == 0:
            print(f"No raw matrix files found, but command line argument -m set to False")

    @staticmethod
    def int_iced_matrix_files(iced_matrixfiles):
        """Round floats in ICE-normalized matrix files to integers."""
        inted_iced_matrixfiles = []

        # Create output directory
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "inted_matrixfiles")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        # Round floats to integers
        for iced_matrixfile in iced_matrixfiles:
            with open(iced_matrixfile) as f:
                lines = f.readlines()

                new_lines = []
                for line in lines:
                    cols = line.strip().split()
                    cols[2] = str(round(float(cols[2])))
                    new_line = '\t'.join(cols) + "\n"
                    new_lines.append(new_line)

                file_name = os.path.basename(iced_matrixfile)

                # Create the output file path with the same name as the original file
                output_file_path = os.path.join(output_dir, file_name)

                # Save the modified matrix file to the new directory
                with open(output_file_path, "w") as f_out:
                    f_out.writelines(new_lines)

                # Add the output file path to the inted_iced_matrixfiles list
                inted_iced_matrixfiles.append(output_file_path)

        return inted_iced_matrixfiles



    @staticmethod
    def group_files_by_experiment_and_resolution(matrix_files, bedfiles):
        grouped_files = {}
        for matrixfile in matrix_files:
            experiment, resolution = FileNameFormatter.extract_experiment_and_resolution(matrixfile)
            key = f"{experiment, resolution}"
            matched_bed_files = [bedfile for bedfile in bedfiles if FileNameFormatter.match_bed_and_matrix_file(bedfile, matrixfile)]
            grouped_files[key] = (matched_bed_files, matrixfile)
        return grouped_files

    @staticmethod
    def group_files(*dirs):
        bedfiles, matrixfiles, iced_matrixfiles = FileManager.find_files(*dirs)

        if SetDirectories.get_normalized_data():
            inted_iced_matrixfiles = FileManager.int_iced_matrix_files(iced_matrixfiles)
            return FileManager.group_files_by_experiment_and_resolution(inted_iced_matrixfiles, bedfiles)
        else:
            return FileManager.group_files_by_experiment_and_resolution(matrixfiles, bedfiles)


# class PipelineInput:
#
#     def __init__(self, root_directories: type(pl.Path), raw_directory: type(str) = "raw", iced_directory: type(str) = "iced", resolutions: type(list[int]) = None):
#         self.root_directories = root_directories
#         self.raw_directory = raw_directory
#         self.iced_directory = iced_directory
#         self.resolutions = resolutions

class PipelineInput:

    def __init__(self, config_loader: ConfigLoader, resolutions: List[int] = None):
        self.input_directories: pl.Path = pl.Path(config_loader.get_value_by_keys("paths", "input_dir"))
        self.output_directory: pl.Path = pl.Path(config_loader.get_value_by_keys("paths", "output_dir"))
        self.raw_directory: str = config_loader.get_value_by_keys("settings", "hicpro_raw_dirname")
        self.iced_directory: str = config_loader.get_value_by_keys("settings", "hicpro_iced_dirname")
        self.tmp_directory: pl.Path = pl.Path(config_loader.get_value_by_keys("paths", "tmp_dir"))
        self.resolutions: List[int] = resolutions
        self.temp_directory = self._get_or_create_tmp_directory(self.output_directory)

        
    def _get_or_create_tmp_directory(self, output_dir: pl.Path) -> pl.Path:
        tmp_dir = output_dir / "tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)  # creates tmp if it doesn't exist
        return tmp_dir

    def _find_files(self) -> Tuple[List[pl.Path], List[pl.Path], List[pl.Path]]:
        bedfiles: List[pl.Path] = []
        matrixfiles: List[pl.Path] = []
        iced_matrixfiles: List[pl.Path] = []

        # Helper function for filtering files
        def filter_files_on_resolution(input_files: List[pl.Path], found_resolutions_in: Set[int]) -> Tuple[List[pl.Path], Set[int]]:
            filtered_files: List[pl.Path] = []
            for file_path in input_files:
                resolution_match = re.search(r'_(\d+)[_.]', file_path.name)
                if resolution_match:
                    resolution = int(resolution_match.group(1))
                    if self.resolutions is None or resolution in self.resolutions:
                        filtered_files.append(file_path)
                    found_resolutions_in.add(resolution)
            return filtered_files, found_resolutions_in

        found_resolutions: Set[int] = set()

        raw_subdirectories: List[pl.Path] = [path for path in self.input_directories.rglob(self.raw_directory) if path.is_dir()]
        iced_subdirectories: List[pl.Path] = [path for path in self.input_directories.rglob(self.iced_directory) if path.is_dir()]

        # Search for bed and matrix files in raw subdirectories
        for subdirectory_path in raw_subdirectories:
            bed_files = list(subdirectory_path.glob('*.bed'))
            matrix_files = list(subdirectory_path.glob('*.matrix'))
            filtered_bed_files, found_resolutions = filter_files_on_resolution(bed_files, found_resolutions)
            bedfiles.extend(filtered_bed_files)
            filtered_matrix_files, found_resolutions = filter_files_on_resolution(matrix_files, found_resolutions)
            matrixfiles.extend(filtered_matrix_files)

        # Search for matrix files in iced subdirectories
        for subdirectory_path in iced_subdirectories:
            iced_matrix_files = list(subdirectory_path.glob('*.matrix'))
            filtered_iced_matrix_files, found_resolutions = filter_files_on_resolution(iced_matrix_files, found_resolutions)
            iced_matrixfiles.extend(filtered_iced_matrix_files)

        file_counts = {
            'bedfiles': len(bedfiles),
            'matrixfiles': len(matrixfiles),
            'iced_matrixfiles': len(iced_matrixfiles)
        }

        if all(count == 0 for count in file_counts.values()) and self.resolutions is not None:
            print(f"No files found for the provided resolutions: {self.resolutions}. "
                  f"Resolutions found in the data: {found_resolutions}")

        # If these conditions are based on certain command-line arguments,
        # make sure to adjust them accordingly in your new structure
        if file_counts["iced_matrixfiles"] == 0:
            print(f"No ICE-normalized matrix files found")
        elif file_counts["matrixfiles"] == 0:
            print(f"No raw matrix files found")

        return bedfiles, matrixfiles, iced_matrixfiles

    def _round_floats_in_iced_files(self, iced_matrixfiles: List[pl.Path]) -> List[pl.Path]:
        inted_iced_matrixfiles = []

        output_dir = self.temp_directory / "inted_matrixfiles"
        output_dir.mkdir(exist_ok=True)

        for iced_matrixfile in iced_matrixfiles:
            with open(iced_matrixfile, 'r') as f:
                new_lines = [
                    "\t".join([cols[0], cols[1], str(round(float(cols[2])))]) + "\n"
                    for cols in (line.strip().split() for line in f)
                ]

            output_file_path = output_dir / iced_matrixfile.name

            with open(output_file_path, "w") as f_out:
                f_out.writelines(new_lines)

            inted_iced_matrixfiles.append(output_file_path)

        return inted_iced_matrixfiles

    def group_files(self) -> Dict[str, Tuple[pl.Path, pl.Path]]:
        bedfiles, matrixfiles, iced_matrixfiles = self._find_files()

        if config_loader.get_value_by_keys("settings", "norm") and iced_matrixfiles:
            inted_iced_matrixfiles = self._round_floats_in_iced_files(iced_matrixfiles)
        else:
            inted_iced_matrixfiles = []

        grouped_raw_files = {}
        grouped_iced_files = {}

        # Extract resolution and experiment name from raw file path
        for matrixfile in matrixfiles:
            resolution = matrixfile.parent.name
            experiment = matrixfile.parents[2].name
            key = f"{experiment, resolution}"

            # Group raw bed file to raw matrix file
            for bedfile in bedfiles:
                if bedfile.stem == matrixfile.stem:
                    grouped_raw_files[key] = (bedfile, matrixfile)

        # Extract resolution and experiment name from ICE-normalized file path
        for inted_matrixfile in inted_iced_matrixfiles:
            file_name = inted_matrixfile.stem
            experiment, resolution, _ = file_name.rsplit("_", 2)
            resolution = int(resolution)
            key = f"{experiment, resolution}"

            # Group ICE-normalized matrix file
            for bedfile in bedfiles:
                bedfile_name = bedfile.stem
                bedfile_experiment, bedfile_resolution, _ = bedfile_name.rsplit("_", 2)
                bedfile_resolution = int(bedfile_resolution)

                if bedfile_experiment == experiment and bedfile_resolution == resolution:
                    grouped_iced_files[key] = (bedfile, inted_matrixfile)

        return grouped_iced_files if SetDirectories.get_normalized_data() else grouped_raw_files

