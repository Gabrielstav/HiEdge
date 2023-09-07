# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules


class blacklist_filter:

    @staticmethod
    def remove_blacklisted_regions(bedpe_file):
        """
        Removes blacklisted regions from bedpe file
        """
        try:
            reference_dir = SetDirectories.get_reference_dir()
            blacklisted_regions_path = os.path.join(reference_dir, "hg19/hg19-blacklist.v2.bed")
            with open(blacklisted_regions_path, "r") as f:
                blacklisted_regions = f.readlines()

            # Filter out mitochondrial DNA and Y chromosome
            with open(bedpe_file, "r") as f:
                bedpe_data = f.readlines()
            filtered_bedpe_data = [line for line in bedpe_data if not any(chrom in line for chrom in ["chrM", "chrY"])]

            blacklisted_pbt = pbt.BedTool(blacklisted_regions)
            blacklisted_bedpe = pbt.BedTool(filtered_bedpe_data)

            window_size = int(re.search(r"(\d+)[^/\d]*$", bedpe_file).group(1))
            if not isinstance(window_size, int):
                raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

            no_overlap_bedpe = blacklisted_bedpe.window(blacklisted_pbt, w=int(window_size), r=False, v=True)
            custom_print(f"Finished processing file: {os.path.basename(bedpe_file)}, PID: {os.getpid()}, TID: {threading.get_ident()}")

            return no_overlap_bedpe

        except Exception as e:
            tid = threading.get_ident()
            print(f"\nError processing file: {bedpe_file}, {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_remove_blacklist():
        """
        Calls the remove_blacklisted_regions function on each file in the BEDPE directory
        """

        global first_print
        first_print = True

        bedpe_dir = SetDirectories.get_temp_dir() + "/bedpe"
        bedpe_files = [os.path.join(bedpe_dir, file) for file in os.listdir(bedpe_dir)]

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "blacklisted")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        try:
            pickle.dumps(bedpe_files)
        except Exception as e:
            print(f"Pickling error: {e}")

        for file in bedpe_files:
            if not os.path.isfile(file):
                print(f"File {file} does not exist.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            futures = list(executor.map(Pipeline.remove_blacklisted_regions, bedpe_files))
            for bedpe_file, future in zip(bedpe_files, futures):
                try:
                    no_blacklist_bedpe = future
                    output_filename = os.path.basename(bedpe_file)[:-len('.bedpe')] + '_no_blacklist.bedpe'
                    converted_nb_bedpe = no_blacklist_bedpe.to_dataframe().to_csv(sep="\t", index=False, header=False)
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        f.writelines(converted_nb_bedpe)

                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")