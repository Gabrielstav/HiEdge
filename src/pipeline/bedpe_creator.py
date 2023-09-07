# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
import pathlib as pl

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