# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from pathlib import Path
from typing import Dict, Set, Tuple
from src.setup.config_loader import ConfigLoader


# Use numpy and pandas, and also validation format (abstract validation, or use dataclass?)
# Then validate data between each step, and do transformtions on dataframes or matrices
# This means might not have to save to file, might be faster, but too memeory intensive, so look into using generators to process data in chunks
# Then implement the stat model using numpy and pandas, and then use generators to process data in chunks as well,
# and look into how this integrates with Cpython, multiprocessing and the GIL
#

class parse_paths:
    pass

class file_formatting:
    pass

class make_bedpe:
    pass

class validate_bedpe:
    pass



class bedpe_creator:

    # Just use dataclass instead ?

    def __init__(self, config_loader: ConfigLoader):
        self.config_loader = config_loader
        self.grouped_files: Dict[str, Tuple[Path, Path]] = {} # Grouped on resolution as key and tuple of bed and matrix files as value
        self.output_directory: Path = Path(config_loader.get_value_by_keys("paths", "output_dir"))
        self.tmp_directory: Path = Path(config_loader.get_value_by_keys("paths", "tmp_dir"))

    def make_bedpe(self, grouped_files: Dict[str, Tuple[Path, Path]]):
        """
        Makes bedpe file from HiC-Pro output
        :param bed_file: BED file from HiC-Pro output grouped by src.setup.pipeline_input.group_files()
        :param matrix_file: Matrix file from HiC-Pro output grouped by src.setup.pipeline_input.group_files()
        :return: BEDPE file saved to temp directory, containing matrix and BED data
        """
        # Main public method to make BEDPE files
        bedpe_creator._parse_paths(self.grouped_files)
        bedpe_creator._file_formatting(self.grouped_files)
        bedpe_creator._make_bedpe(self.grouped_files)

    def _parse_paths(self):



    def _file_formatting(self):
        pass


class BEDPE_Creator:

    @staticmethod
    def make_bedpe(bed_file, matrix_file):
        """
        Makes bedpe file from HiC-Pro output
        """

        try:
            bedpe = []

            bed_lines = open(bed_file, "r").readlines()
            matrix_lines = open(matrix_file, "r").readlines()

            bed_dict = {}
            for line in bed_lines:
                line = line.strip().split("\t")
                bed_dict[line[3]] = line
            for line in matrix_lines:
                line = line.strip().split("\t")
                bedpe.append(f"{bed_dict[line[0]][0]}"
                             f"\t{bed_dict[line[0]][1]}"
                             f"\t{bed_dict[line[0]][2]}"
                             f"\t{bed_dict[line[1]][0]}"
                             f"\t{bed_dict[line[1]][1]}"
                             f"\t{bed_dict[line[1]][2]}"
                             f"\t{line[2]}\n")
            return bedpe

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error processing file: {bed_file, matrix_file}, {e}, PID: {os.getpid()}, TID: {tid}")
            raise


    @staticmethod
    def input_to_make_bedpe(grouped_files):
        """
        Makes bedpe files from HiC-Pro output and saves them to the temp directory
        :param grouped_files: Output of Pipeline_Input.group_files()
        :return: BEDPE files saved to temp directory
        """

        global first_print
        first_print = True

        os.chdir(SetDirectories.get_temp_dir())
        if not os.path.exists("bedpe"):
            os.mkdir("bedpe")
        else:
            shutil.rmtree("bedpe")
            os.mkdir("bedpe")
        os.chdir("bedpe")

        for key, val in grouped_files.items():

            split_key = key.split(",")
            experiment = split_key[0].replace("'", "").replace("(", "").replace(")", "").replace(" ", "_")
            resolution = split_key[1].replace("'", "").replace("(", "").replace(")", "").replace(" ", "_")
            # Remove trailing underscore or space
            if experiment[-1] in ("_", " "):
                experiment = experiment[:-1]
            if resolution[-1] in ("_", " "):
                resolution = resolution[:-1]
            # Replace consecutive underscores with a single underscore
            experiment = '_'.join(filter(None, experiment.split('_')))
            resolution = '_'.join(filter(None, resolution.split('_')))

            bedfile = val[0]
            matrixfile = val[1]
            bedpe = Pipeline.make_bedpe(bedfile, matrixfile)

            with open(experiment + "_" + resolution + ".bedpe", "w") as f:
                f.writelines(bedpe)
                f.close()