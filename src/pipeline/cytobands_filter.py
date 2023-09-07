# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules

class cytobands_filter:

    @staticmethod
    def remove_cytobands(blacklisted_bedpe_file):
        """
        Cytoband locations are determined in this case by Giemsa staining
        and are located and removed from the BEDPE file.
        """

        try:
            os.chdir(SetDirectories.get_reference_dir())
            cytobands = os.path.join(SetDirectories.get_reference_dir(), "hg19/cytoBand_hg19.txt")
            centromeric_regions = []
            with open(cytobands, "r") as f:
                cytobands = f.readlines()
                for line in cytobands:
                    line = line.strip().split("\t")
                    if line[4] == "acen":
                        centromeric_regions.append(line[0:5])

            os.chdir(SetDirectories.get_temp_dir() + "/blacklisted")
            blacklisted_pbt = pbt.BedTool(blacklisted_bedpe_file)

            window_size = int(re.search(r"(\d+)[^/\d]*$", blacklisted_bedpe_file).group(1))
            if not isinstance(window_size, int):
                raise ValueError(f"Window size must be an integer, {window_size} is not an integer.")

            no_cytobands = blacklisted_pbt.window(centromeric_regions, w=int(window_size), r=False, v=True)
            custom_print(f"Finished processing file: {os.path.basename(blacklisted_bedpe_file)}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return no_cytobands

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error processing {blacklisted_bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_remove_cytobands():
        """
        Calls the remove_cytobands function on each file in the blacklisted directory
        """

        global first_print
        first_print = True

        blacklisted_dir_path = SetDirectories.get_temp_dir() + "/blacklisted"
        blacklisted_dir = os.listdir(blacklisted_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "no_cytobands")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for file in blacklisted_dir:
            full_path = os.path.join(blacklisted_dir_path, file)
            if not os.path.isfile(full_path):
                print(f"File {file} does not exist.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            full_paths = [os.path.join(blacklisted_dir_path, file) for file in blacklisted_dir]
            futures = list(executor.map(Pipeline.remove_cytobands, full_paths))
            for bedpe_file, future in zip(blacklisted_dir, futures):
                try:
                    no_cytoband_bedpe = future
                    output_filename = f"{bedpe_file[:-len('.bedpe')]}_no_cytobands.bedpe"
                    converted_nc_bedpe = no_cytoband_bedpe.to_dataframe().to_csv(sep="\t", index=False, header=False)
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        f.writelines(converted_nc_bedpe)

                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {bedpe_file}: {e}, PID: {os.getpid()}, TID: {tid}")