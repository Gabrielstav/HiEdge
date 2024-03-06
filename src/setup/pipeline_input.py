# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import re
from pathlib import Path
from typing import List, Set, Tuple, Iterable, Optional
from src.setup.config_loader import Config
from src.setup.containers import GroupedFiles, Metadata

class HicProInputFilePreparer:

    def __init__(self, config: Config):
        self.config = config
        self.file_finder = FileFinder(config)
        self.file_grouper = HicProFileGrouper(config)
        # self.matrix_rounder = IcedMatrixRounder(config)

    def prepare_input_files(self) -> List[GroupedFiles]:
        bedfiles, matrixfiles, biasfiles = self.file_finder.find_hicpro_bed_matrix_files()
        grouped_files = self.file_grouper.group_files(bedfiles, matrixfiles, biasfiles)
        return grouped_files

class FileFinder:

    def __init__(self, config: Config):
        self.config = config
        if config.pipeline_settings.interaction_type == "intra":
            self.desired_resolutions = set(config.pipeline_settings.intra_resolutions)
        elif config.pipeline_settings.interaction_type == "inter":
            self.desired_resolutions = set(config.pipeline_settings.inter_resolutions)
        elif config.pipeline_settings.interaction_type == "mixed":
            self.desired_resolutions = set(config.pipeline_settings.intra_resolutions + config.pipeline_settings.inter_resolutions)
        else:
            raise ValueError(f"Unknown interaction type: {config.pipeline_settings.interaction_type}")

    def find_hicpro_bed_matrix_files(self) -> Tuple[List[Path], List[Path], Optional[List[Path]]]:
        """
        Discovers bed files and matrix files (either raw or iced) based on the normalization setting.
        :returns: Tuple of lists of Paths to the bed files and matrix files.
        """
        bed_files: List[Path] = []
        selected_matrix_files: List[Path] = []
        found_resolutions: Set[int] = set()

        raw_subdirectories = self._find_subdirectories(self.config.pipeline_settings.hicpro_raw_dirname)
        iced_subdirectories = self._find_subdirectories(self.config.pipeline_settings.hicpro_norm_dirname)

        # TODO: remove later
        print(f"Searching in raw subdirectories: {raw_subdirectories}")
        print(f"Found bed files: {bed_files}")
        print(f"Found matrix files: {selected_matrix_files}")


        # Search for bed and matrix files in raw subdirectories
        for subdir in raw_subdirectories:
            found_bed_files, found_resolutions = self._filter_files_on_resolution(subdir.glob("*.bed"), found_resolutions)
            bed_files.extend(found_bed_files)

            if not self.config.pipeline_settings.iced_data:
                found_matrix_files, found_resolutions = self._filter_files_on_resolution(subdir.glob("*.matrix"), found_resolutions)
                selected_matrix_files.extend(found_matrix_files)

        # If iced_matrix config, only search in iced subdirectories
        if self.config.pipeline_settings.iced_data:
            for subdir in iced_subdirectories:
                found_iced_matrix_files, found_resolutions = self._filter_files_on_resolution(subdir.glob("*.matrix"), found_resolutions)
                selected_matrix_files.extend(found_iced_matrix_files)

        bias_files = []
        if self.config.statistical_settings.use_hicpro_bias:
            for subdir in iced_subdirectories:
                found_files, _ = self._filter_files_on_resolution(subdir.glob("*.biases"))
                bias_files.extend(found_files)

        # Check for missing resolutions after finding files
        missing_resolutions = self.desired_resolutions - found_resolutions
        for resolution in missing_resolutions:
            print(f"Data for resolution {resolution} not present in input directory.")

        if missing_resolutions:
            for resolution in missing_resolutions:
                print(f"Data for resolution {resolution} not present in input directory.")
        else:
            print("Successfully found input files for all specified resolutions.")

        return bed_files, selected_matrix_files, bias_files

    def _find_subdirectories(self, dirname: str) -> List[Path]:
        return [path for path in self.config.paths.input_dir.glob(f"**/{dirname}/*") if path.is_dir()]

    def _filter_files_on_resolution(self, input_files: Iterable[Path], found_resolutions=None) -> Tuple[List[Path], Set[int]]:
        """
        Filters input files based on their resolution.

        :returns: Tuple of filtered Paths and the set of found resolutions.
        """
        if found_resolutions is None:
            found_resolutions = set()
        filtered_files: List[Path] = []
        for file_path in input_files:
            resolution_match = re.search(r"_(\d+)[_.]", file_path.name)
            if resolution_match:
                resolution = int(resolution_match.group(1))
                print(f"Found resolution: {resolution} for file: {file_path}")
                if self.desired_resolutions is None or resolution in self.desired_resolutions:
                    filtered_files.append(file_path)
                found_resolutions.add(resolution)
        return filtered_files, found_resolutions


class HicProFileGrouper:

    def __init__(self, config: Config):
        self.config = config
        self.desired_resolutions = set(config.pipeline_settings.intra_resolutions + config.pipeline_settings.inter_resolutions)

    def group_files(self, bedfiles: List[Path], matrixfiles: List[Path], biasfiles: Optional[List[Path]] = None) -> List[GroupedFiles]:

        """
        Groups bed files and matrix files based on the normalization setting.

        :param bedfiles: List of Paths to the bed files.
        :param matrixfiles: List of Paths to the matrix files.
        :param biasfiles: List of Paths to the bias files.
        :returns: Dictionary with keys of the form (experiment, resolution) and values of the form (bedfile, matrixfile).
        """

        def extract_metadata_from_file(file: Path) -> Tuple[str, int]:
            """
            Extracts the experiment and resolution from a given file based on its name and path.
            """
            if self.config.pipeline_settings.iced_data:
                exp_name, res_value, _ = file.stem.rsplit("_", 2)
                res_value = int(res_value)
            else:
                exp_name = file.parents[2].name
                res_value = int(file.stem.rsplit("_", 2)[1])
            return exp_name, res_value

        grouped_files = []

        bedfile_lookup = {extract_metadata_from_file(bedfile): bedfile for bedfile in bedfiles}
        biasfile_lookup = {extract_metadata_from_file(biasfile): biasfile for biasfile in biasfiles or []}

        for matrixfile in matrixfiles:
            key = extract_metadata_from_file(matrixfile)
            if key in bedfile_lookup:
                experiment, resolution = key
                bias_file = biasfile_lookup.get(key) if self.config.statistical_settings.use_hicpro_bias else None
                metadata = Metadata(experiment=experiment, resolution=resolution, bias_file_path=bias_file)
                grouped_files.append(GroupedFiles(metadata=metadata, bed_file=bedfile_lookup[key], matrix_file=matrixfile))

        for resolution in self.desired_resolutions:
            groups_for_resolution = [group for group in grouped_files if group.metadata.resolution == resolution]
            if not groups_for_resolution:
                print(f"No files found for resolution {resolution}.")
            else:
                if self.config.statistical_settings.use_hicpro_bias and any(group.metadata.bias_file_path is None for group in groups_for_resolution):
                    print(f"Bias files missing for resolution {resolution}.")

        return grouped_files
