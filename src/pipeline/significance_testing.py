# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules
from scipy import stats

# stats.nchypergeom_wallenius
# stats.nchypergeom_fisher
# stats.nbinom
from typing import List




class NchgDistribution:

    def __init__(self):
        pass

class Splines:

    def __init__(self):
        pass

class AggregateStatistics:

    """
    this class calculates aggregate statistics over every chunk (total interaction count, expected interaction count)
    """

    def __init__(self):
        pass

class FindSignificantInteractions:

    """
    this class should use the NCHG distribution and the splines, as well as aggregate statistics for each chunk (by chromosome?),
    and apply the significance testing to each chunk (each interaction needs genomic distance, observed interaction count + all aggregate statistics)
    """

    def __init__(self):
        pass

class FindSignificantInterchromosomalInteractions:

    """
    For interchromosomal interactions we do not needd splines nor the NCHG distribution, because there is no distance dependant decay,
    we use NCHG dist and Fiscer's exact test (?). This will be more difficult to parallelize, but can chunk by source interaction
    to X amount of targets across chromosomes, so each chunks contains all interactions for one chromosome, with summary statistics derived from same method as the intra data.
    """

    def __init__(self):
        pass


class significant_contacts:

    @staticmethod
    def find_significant_interactions(bedpe_file):
        """
        NCHG script to calculate the significance of interactions:
        m = minimum interaction length in bp, should be same as window size used to make the bedpe file (resolution)
        p = print counts, needed for FDR correction for padj with BH method.
        i = Use interchromosomal interactions as well when calculating p-values (w in args for this script)
        """

        try:
            # Setting min interactions length to same as bin size from HiC-Pro
            window_size = int(re.search(r"_(\d+)_", bedpe_file).group(1))  # Any number in file name with underscores on both sides
            if not isinstance(window_size, int):
                raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

            # Run NCHG
            nchg_flags = []

            if SetDirectories.inter_interactions:
                nchg_flags.append("-i")

            elif SetDirectories.mixed_interactions:
                intrachromosomal = True
                with open(bedpe_file, "r") as f:
                    for line in f:
                        columns = line.strip().split("\t")
                        if columns[0] != columns[3]:
                            intrachromosomal = False
                            break

                if not intrachromosomal:
                    nchg_flags.append("-i")

            if SetDirectories.n_quantiles:
                nchg_flags.append("-n ")
                nchg_flags.append(str(SetDirectories.n_quantiles))

            if not SetDirectories.inter_interactions and not SetDirectories.intra_interactions and not SetDirectories.mixed_interactions:
                raise ValueError("Must select at least one type of interaction type for significance testing: Inter, intra, or mixed.")

            # Run NCHG and clean up completed processes
            nchg_command = [SetDirectories.get_nchg_path(), bedpe_file, "-m", str(window_size), "-p"] + nchg_flags
            nchg_run = sp.Popen(nchg_command, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = nchg_run.communicate()
            nchg_run.kill()
            nchg_run.wait()

            print(f"Finished processing file: {os.path.basename(bedpe_file)}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return stdout.decode("utf-8").split("\t")

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def split_bedpe_by_chromosome(bedpe_file, output_dir):
        """
        For splitting intrachromosomal input files to NCHG, gives different results than not splitting, but much faster.
        Split a bedpe file by chromosome, only used for intra data split on chromosome.
        Does not currently work for inter data, since that needs the whole genome for statistical significance testing.
        :param bedpe_file: BEDPE file with interactions and blacklisted/centromeric regions removed
        :param output_dir: Directory containing the split BEDPE files, one BEDPE file per chromosome
        :return: BEDPE file containing one chromosome's worth of interactions
        """
        try:
            with open(bedpe_file, "r") as f:
                data = f.readlines()

            chromosomes = {}
            for line in data:
                chr_name = line.split("\t")[0]
                if chr_name not in chromosomes:
                    chromosomes[chr_name] = []
                chromosomes[chr_name].append(line.strip())  # Strip whitespace from the line

            # Write each chromosome to a separate file
            chr_files = []
            for chr_name, chr_data in chromosomes.items():
                output_filename = f"{os.path.basename(bedpe_file)[:-len('.bedpe')]}_{chr_name}_split.bedpe"
                output_filepath = os.path.join(output_dir, output_filename)
                with open(output_filepath, "w") as f:
                    f.write("\n".join(chr_data) + "\n")
                chr_files.append(output_filepath)

            return chr_files

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise



    @staticmethod
    def split_bedpe_by_chromosome_pairs(bedpe_file, output_dir):
        """
        For splitting input files to NCHG in interchromosomal interactions. Gives different results than not splitting, but much faster.
        Split a bedpe file by chromosome pairs, only used for inter data split on chromosome pairs.
        This gives inter interactions, but needs to be validated against no_split inter interactions.
        :param bedpe_file: BEDPE file containing inter interactions with blacklisted/centromeric regions removed.
        :param output_dir: Output directory for the split BEDPE files, one BEDPE file per chromosome pair.
        :return: BEDPE file containing one chromosome pair's worth of inter interactions.
        """
        try:
            with open(bedpe_file, "r") as f:
                data = f.readlines()

            # Group data by chromosome pairs
            chromosome_pairs = {}
            for line in data:
                columns = line.split("\t")
                chr1, chr2 = columns[0], columns[3]
                chr_pair_key = tuple(sorted([chr1, chr2]))

                if chr_pair_key in chromosome_pairs:
                    chromosome_pairs[chr_pair_key].append(line.strip())
                else:
                    chromosome_pairs[chr_pair_key] = [line.strip()]

            # Write each chromosome pair to a separate file
            chr_files = []
            for (chr1, chr2), chr_data in chromosome_pairs.items():
                output_filename = f"{os.path.basename(bedpe_file)[:-len('.bedpe')]}_{chr1}-{chr2}_split.bedpe"
                output_filepath = os.path.join(output_dir, output_filename)
                with open(output_filepath, "w") as f:
                    f.write("\n".join(chr_data) + "\n")
                chr_files.append(output_filepath)

            return chr_files

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_nchg():
        """
        Calls the NCHG script on all files in the no_cytobands directory to find significant interactions
        """

        no_cytobands_dir_path = SetDirectories.get_temp_dir() + "/no_cytobands"
        no_cytobands_dir = os.listdir(no_cytobands_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "NCHG_output")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        if SetDirectories.get_no_split():
            all_files = [os.path.join(no_cytobands_dir_path, file) for file in no_cytobands_dir]
            input_file_map = {file: file for file in all_files}
            files_to_process = all_files

        elif SetDirectories.inter_interactions:
            # Create a directory to store split chromosome files
            chr_split_base_dir = os.path.join(SetDirectories.get_temp_dir(), "input_to_nchg")
            if not os.path.exists(chr_split_base_dir):
                os.mkdir(chr_split_base_dir)
            else:
                shutil.rmtree(chr_split_base_dir)
                os.mkdir(chr_split_base_dir)

            # Split input files by chromosome pairs
            all_chr_files = []
            # keep original input file name for each chromosome file
            input_file_map = {}

            for file in no_cytobands_dir:
                full_path = os.path.join(no_cytobands_dir_path, file)
                if not os.path.isfile(full_path):
                    print(f"File {file} does not exist.")

                # Create a subdirectory for each input file inside the chr_split_base_dir
                file_split_dir = os.path.join(chr_split_base_dir, f"{file[:-len('_no_blacklist_no_cytobands.bedpe')]}_split")
                if not os.path.exists(file_split_dir):
                    os.mkdir(file_split_dir)
                else:
                    shutil.rmtree(file_split_dir)
                    os.mkdir(file_split_dir)

                chr_files = Pipeline.split_bedpe_by_chromosome_pairs(full_path, file_split_dir)

                for chr_file in chr_files:
                    input_file_map[chr_file] = file
                all_chr_files.extend(chr_files)

            files_to_process = all_chr_files

        else:
            # Create a directory to store split chromosome files
            chr_split_base_dir = os.path.join(SetDirectories.get_temp_dir(), "input_to_nchg")
            if not os.path.exists(chr_split_base_dir):
                os.mkdir(chr_split_base_dir)
            else:
                shutil.rmtree(chr_split_base_dir)
                os.mkdir(chr_split_base_dir)

            # Split input files by chromosome
            all_chr_files = []
            # keep original input file name for each chromosome file
            input_file_map = {}

            for file in no_cytobands_dir:
                full_path = os.path.join(no_cytobands_dir_path, file)
                if not os.path.isfile(full_path):
                    print(f"File {file} does not exist.")

                # Create a subdirectory for each input file inside the chr_split_base_dir
                file_split_dir = os.path.join(chr_split_base_dir, f"{file[:-len('_no_blacklist_no_cytobands.bedpe')]}_split")
                if not os.path.exists(file_split_dir):
                    os.mkdir(file_split_dir)
                else:
                    shutil.rmtree(file_split_dir)
                    os.mkdir(file_split_dir)

                chr_files = Pipeline.split_bedpe_by_chromosome(full_path, file_split_dir)

                for chr_file in chr_files:
                    input_file_map[chr_file] = file
                all_chr_files.extend(chr_files)

            files_to_process = all_chr_files

        # Run find_significant_interactions on chromosome-specific files in parallel
        output_file_data = defaultdict(list)
        executorclass = concurrent.futures.ProcessPoolExecutor if executor_type == 'multiprocessing' else concurrent.futures.ThreadPoolExecutor
        with executorclass(max_workers=SetDirectories.get_threads()) as executor:
            futures = list(executor.map(Pipeline.find_significant_interactions, files_to_process))
            for bedpe_file, future in zip(files_to_process, futures):
                try:
                    nchg_output = future
                    input_file = input_file_map[bedpe_file]
                    output_file_data[input_file].append(nchg_output)

                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")

        executor.shutdown(wait=True)

        # Merge output files back together
        for input_file, nchg_outputs in output_file_data.items():
            output_filename = f"{os.path.basename(input_file)[:-len('_no_blacklist_no_cytobands.bedpe')]}_nchg_output.txt"
            output_filepath = os.path.join(output_dir, output_filename)
            with open(output_filepath, "w") as f:
                for nchg_output in nchg_outputs:
                    f.writelines(nchg_output)