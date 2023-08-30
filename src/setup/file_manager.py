# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

import os
import re
import shutil as shutil


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
            inted_iced_matrixfiles = FileManager.intify_iced_matrix_files(iced_matrixfiles)
            return FileManager.group_files_by_experiment_and_resolution(inted_iced_matrixfiles, bedfiles)
        else:
            return FileManager.group_files_by_experiment_and_resolution(matrixfiles, bedfiles)
