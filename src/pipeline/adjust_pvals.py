# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules

class adjust_pvals:

    @staticmethod
    def adjust_pvalues(nchg_file, fdr_thresh=SetDirectories.get_fdr_threshold(), log_ratio_threshold=2, method="fdr_bh"):
        """
        Adjusts the p-values using the Benjamini-Hochberg method
        """

        try:
            # Finds the p-values and log ratios of the interactions
            p_values = []
            processed_lines = []
            with open(nchg_file, "r") as nchg_f:
                for line in nchg_f:
                    col = line.split()
                    if len(col) < 11:  # Expect 11 columns in the input file
                        print(f"Error: {nchg_file} does not have the correct number of columns.")
                    if col[7] == "0":
                        continue  # Skips the line if interactions/edges is 0

                    p_values.append(float(col[6]))

                    log_ratio = 0.0
                    if float(col[9]) != 0 and float(col[10]) != 0:
                        log_ratio = math.log(float(col[9]), 2) - math.log(float(col[10]), 2)

                    processed_line = ' '.join(col) + ' ' + str(log_ratio)
                    processed_lines.append(processed_line)

            padj = list(multicomp.multipletests(p_values, method=method))
            padj_out = []

            # Filters the interactions based on the log ratio and FDR thresholds
            for i, processed_line in enumerate(processed_lines):
                col = processed_line.split()
                col.append(str(padj[1][i]))
                if float(col[11]) >= log_ratio_threshold and float(col[12]) <= fdr_thresh:
                    padj_out.append("\t".join(col[:6] + [col[12]]))

            print(f"Finished processing file: {os.path.basename(nchg_file)}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return padj_out

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {nchg_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_adjust_pvalues():
        """
        Calls the adjust_pvalues function on all files in the NCHG_output directory
        """

        nchg_dir_path = SetDirectories.get_temp_dir() + "/NCHG_output"
        nchg_dir = os.listdir(nchg_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "padj")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for file in nchg_dir:
            full_path = os.path.join(nchg_dir_path, file)
            if not os.path.isfile(full_path):
                print(f"File {file} does not exist.")
                continue

        with concurrent.futures.ProcessPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            full_paths = [os.path.join(nchg_dir_path, file) for file in nchg_dir]
            futures = list(executor.map(Pipeline.adjust_pvalues, full_paths))
            for nchg_file, future in zip(nchg_dir, futures):
                try:
                    padj = future
                    output_filename = f"{nchg_file[:-len('_nchg_output.txt')]}_padj.txt"
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        for line in padj:
                            f.write(line + "\n")
                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {nchg_file}: {e}, PID: {os.getpid()}, TID: {tid}")