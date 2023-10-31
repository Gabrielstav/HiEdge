# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import re
from pathlib import Path
from typing import List, Set, Tuple, Dict, Iterable, Optional
from config_loader import Config


class HicProInputFilePreparer:

    def __init__(self, config: Config):
        self.config = config
        self.file_finder = FileFinder(config)
        self.file_grouper = HicProFileGrouper(config)
        self.matrix_rounder = IcedMatrixRounder(config)

    def prepare_input_files(self) -> Dict[Tuple[str, int], Tuple[Path, Path]]:
        bedfiles, matrixfiles = self.file_finder.find_hicpro_bed_matrix_files()
        grouped_files = self.file_grouper.group_files(bedfiles, matrixfiles)

        if self.config.pipeline_settings.normalized_data:
            matrix_files_to_round = [value[1] for value in grouped_files.values()]
            rounded_matrix_files = self.matrix_rounder.round_floats_in_iced_files(matrix_files_to_round)

            # Replace float matrix files in the grouped_files dict with the inted rounded ones
            for key, (bedfile, _) in grouped_files.items():
                grouped_files[key] = (bedfile, rounded_matrix_files.pop(0))

        return grouped_files

class FileFinder:

    def __init__(self, config: Config):
        self.config = config
        self.resolutions: Optional[List[Path]] = None

    def find_hicpro_bed_matrix_files(self) -> Tuple[List[Path], List[Path]]:
        """
        Discovers bed files and matrix files (either raw or iced) based on the normalization setting.
        :returns: Tuple of lists of Paths to the bed files and matrix files.
        """
        bed_files: List[Path] = []
        selected_matrix_files: List[Path] = []
        found_resolutions: Set[int] = set()

        raw_subdirectories = self._find_subdirectories(self.config.pipeline_settings.hicpro_raw_dirname)
        iced_subdirectories = self._find_subdirectories(self.config.pipeline_settings.hicpro_norm_dirname)

        # Search for bed and matrix files in raw subdirectories
        for subdir in raw_subdirectories:
            found_bed_files, found_resolutions = self._filter_files_on_resolution(subdir.glob('*.bed'), found_resolutions)
            bed_files.extend(found_bed_files)

            if not self.config.pipeline_settings.normalized_data:
                found_matrix_files, found_resolutions = self._filter_files_on_resolution(subdir.glob('*.matrix'), found_resolutions)
                selected_matrix_files.extend(found_matrix_files)

        # If normalized_data_flag, only search in iced subdirectories
        if self.config.pipeline_settings.normalized_data:
            for subdir in iced_subdirectories:
                found_iced_matrix_files, found_resolutions = self._filter_files_on_resolution(subdir.glob('*.matrix'), found_resolutions)
                selected_matrix_files.extend(found_iced_matrix_files)

        return bed_files, selected_matrix_files

    def _find_subdirectories(self, dirname: str) -> List[Path]:
        return [path for path in self.config.paths.input_dir.rglob(dirname) if path.is_dir()]

    def _filter_files_on_resolution(self, input_files: Iterable[Path], found_resolutions=None) -> Tuple[List[Path], Set[int]]:
        """
        Filters input files based on their resolution.

        :returns: Tuple of filtered Paths and the set of found resolutions.
        """
        if found_resolutions is None:
            found_resolutions = set()
        filtered_files: List[Path] = []
        for file_path in input_files:
            resolution_match = re.search(r'_(\d+)[_.]', file_path.name)
            if resolution_match:
                resolution = int(resolution_match.group(1))
                if self.resolutions is None or resolution in self.resolutions:
                    filtered_files.append(file_path)
                found_resolutions.add(resolution)
        return filtered_files, found_resolutions


class HicProFileGrouper:

    def __init__(self, config: Config):
        self.config = config

    def group_files(self, bedfiles: List[Path], matrixfiles: List[Path]) -> Dict[Tuple[str, int], Tuple[Path, Path]]:
        """
        Groups bed files and matrix files based on the normalization setting.

        :param bedfiles: List of Paths to the bed files.
        :param matrixfiles: List of Paths to the matrix files.
        :returns: Dictionary with keys of the form (experiment, resolution) and values of the form (bedfile, matrixfile).
        """

        def extract_metadata_from_file(file: Path) -> Tuple[str, int]:
            """
            Extracts the experiment and resolution from a given file based on its name and path.
            """
            if self.config.pipeline_settings.normalized_data:
                experiment, resolution, _ = file.stem.rsplit("_", 2)
                resolution = int(resolution)
            else:
                experiment = file.parents[2].name
                resolution = int(file.stem.rsplit("_", 2)[1])
            return experiment, resolution

        grouped_files: Dict[Tuple[str, int], Tuple[Path, Path]] = {}

        bedfile_lookup = {extract_metadata_from_file(bedfile): bedfile for bedfile in bedfiles}

        for matrixfile in matrixfiles:
            key = extract_metadata_from_file(matrixfile)
            if key in bedfile_lookup:
                grouped_files[key] = (bedfile_lookup[key], matrixfile)

        return grouped_files

class IcedMatrixRounder:

    def __init__(self, config: Config):
        self.config = config

    def round_floats_in_iced_files(self, iced_matrixfiles: List[Path]) -> List[Path]:
        """
        Rounds the float values in ICE-normalized matrix files and saves them in a new directory
        inside the output directory. The new files are returned.

        :param: List of Paths to the iced matrix files.
        :returns: List of Paths to the rounded iced matrix files.
        """
        inted_iced_matrixfiles = []

        for iced_matrixfile in iced_matrixfiles:
            output_file_path = self.config.paths.output_dir / iced_matrixfile.name

            if output_file_path.exists():
                # Handle existing file, e.g., skip, overwrite, or rename
                continue

            with open(iced_matrixfile, 'r') as f_in, open(output_file_path, "w") as f_out:
                for line in f_in:
                    cols = line.strip().split()
                    rounded_value = round(float(cols[2]))
                    f_out.write("\t".join([cols[0], cols[1], str(rounded_value)]) + "\n")

            inted_iced_matrixfiles.append(output_file_path)

        return inted_iced_matrixfiles
