# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import re
from pathlib import Path
from typing import List, Set, Tuple, Dict
from config_loader import ConfigLoader
from dataclasses import dataclass

@dataclass
class configuration:
    hicpro_input_directory: Path
    output_directory: Path
    raw_directory_name: str
    iced_directory_name: str
    tmp_directory: Path
    resolutions: List[int]
    normalized_data_flag: bool

class file_finder:
    pass

class file_grouper:
    pass

class iced_rounder_intifier:
    pass



class PipelineInputFromHiCPro:

    # Refactor to use dataclass instead ?
    # This would be better because we can use this as type hinting for the return type of the prepare_input_files method instead of the Dict[str, Tuple[Path, Path]] ?

    def __init__(self, config_loader: ConfigLoader, resolutions: List[int] = None):
        self.config_loader = config_loader
        self.input_directory: Path = Path(config_loader.get_value_by_keys("paths", "input_dir"))
        self.output_directory: Path = Path(config_loader.get_value_by_keys("paths", "output_dir"))
        self.raw_directory_name: str = config_loader.get_value_by_keys("settings", "hicpro_raw_dirname")
        self.iced_directory_name: str = config_loader.get_value_by_keys("settings", "hicpro_iced_dirname")
        self.tmp_directory: Path = Path(config_loader.get_value_by_keys("paths", "tmp_dir"))
        self.resolutions: List[int] = resolutions
        self.normalized_data_flag: bool = config_loader.get_value_by_keys("settings", "normalized_data")

    def prepare_input_files(self):
        """
        Main public method to prepare input files for the pipeline.
        """
        bedfiles, matrixfiles = self._find_files()
        grouped_files = self._group_files(bedfiles, matrixfiles)

        if self.normalized_data_flag:
            output_dir = self.tmp_directory / "inted_matrixfiles"
            output_dir.mkdir(exist_ok=True)

            matrix_files_to_round = [value[1] for value in grouped_files.values()]
            rounded_matrix_files = self._round_floats_in_iced_files(matrix_files_to_round)

            # Replace float matrix files in the grouped_files dict with the inted rounded ones
            for key, (bedfile, _) in grouped_files.items():
                grouped_files[key] = (bedfile, rounded_matrix_files.pop(0))

        return grouped_files

    def _find_files(self) -> Tuple[List[Path], List[Path]]:
        """
        Discovers bed files and matrix files (either raw or iced) based on the normalization setting.

        :returns: Tuple of lists of Paths to the bed files and matrix files.
        """
        bedfiles: List[Path] = []
        selected_matrixfiles: List[Path] = []

        def filter_files_on_resolution(input_files: List[Path], found_resolutions_in: Set[int]) -> Tuple[List[Path], Set[int]]:
            filtered_files: List[Path] = []
            for file_path in input_files:
                resolution_match = re.search(r'_(\d+)[_.]', file_path.name)
                if resolution_match:
                    resolution = int(resolution_match.group(1))
                    if self.resolutions is None or resolution in self.resolutions:
                        filtered_files.append(file_path)
                    found_resolutions_in.add(resolution)
            return filtered_files, found_resolutions_in

        found_resolutions: Set[int] = set()
        raw_subdirectories: List[Path] = [path for path in self.input_directory.rglob(self.raw_directory) if path.is_dir()]
        iced_subdirectories: List[Path] = [path for path in self.input_directory.rglob(self.iced_directory) if path.is_dir()]

        # Search for bed and matrix files in raw subdirectories
        for subdirectory_path in raw_subdirectories:
            bed_files = list(subdirectory_path.glob('*.bed'))
            matrix_files = list(subdirectory_path.glob('*.matrix'))
            filtered_bed_files, found_resolutions = filter_files_on_resolution(bed_files, found_resolutions)
            bedfiles.extend(filtered_bed_files)
            if not self.normalized_data_flag:
                filtered_matrix_files, found_resolutions = filter_files_on_resolution(matrix_files, found_resolutions)
                selected_matrixfiles.extend(filtered_matrix_files)

        # If normalized_data_flag, only search in iced subdirectories
        if self.normalized_data_flag:
            for subdirectory_path in iced_subdirectories:
                iced_matrix_files = list(subdirectory_path.glob('*.matrix'))
                filtered_iced_matrix_files, found_resolutions = filter_files_on_resolution(iced_matrix_files, found_resolutions)
                selected_matrixfiles.extend(filtered_iced_matrix_files)

        return bedfiles, selected_matrixfiles

    def _group_files(self, bedfiles: List[Path], matrixfiles: List[Path]) -> Dict[str, Tuple[Path, Path]]:
        """
        Groups bed files and matrix files (either raw or iced) based on the normalization setting.

        :param bedfiles: List of Paths to the bed files.
        :param matrixfiles: List of Paths to the matrix files.
        :returns: Dictionary with keys of the form (experiment, resolution) and values of the form (bedfile, matrixfile).
        """

        grouped_files: Dict[str, Tuple[Path, Path]] = {}

        bedfile_lookup = {}
        for bedfile in bedfiles:
            if self.normalized_data_flag:
                bedfile_experiment, bedfile_resolution, _ = bedfile.stem.rsplit("_", 2)
                bedfile_resolution = int(bedfile_resolution)
                key = (bedfile_experiment, bedfile_resolution)
            else:
                key = bedfile.stem
            bedfile_lookup[key] = bedfile

        for matrixfile in matrixfiles:
            # Check if using iced matrix files or raw matrix files
            resolution = matrixfile.parent.name if not self.normalized_data_flag else matrixfile.stem.rsplit("_", 2)[1]
            experiment = matrixfile.parents[2].name if not self.normalized_data_flag else matrixfile.stem.rsplit("_", 2)[0]

            key = (experiment, int(resolution)) if self.normalized_data_flag else experiment
            if key in bedfile_lookup:
                grouped_files[f"{experiment, resolution}"] = (bedfile_lookup[key], matrixfile)

        return grouped_files

    def _round_floats_in_iced_files(self, iced_matrixfiles: List[Path]) -> List[Path]:
        """
        Rounds the float values in ICE-normalized matrix files and saves them in a new directory
        inside the tmp directory. The new files are returned.

        :param: List of Paths to the iced matrix files.
        :returns: List of Paths to the rounded iced matrix files.
        """
        inted_iced_matrixfiles = []

        for iced_matrixfile in iced_matrixfiles:
            with open(iced_matrixfile, 'r') as f:
                new_lines = [
                    "\t".join([cols[0], cols[1], str(round(float(cols[2])))]) + "\n"
                    for cols in (line.strip().split() for line in f)
                ]

            output_file_path = self.output_directory / iced_matrixfile.name

            with open(output_file_path, "w") as f_out:
                f_out.writelines(new_lines)

            inted_iced_matrixfiles.append(output_file_path)

        return inted_iced_matrixfiles


