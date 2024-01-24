# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

# Import modules

class make_networks:

    @staticmethod
    def make_weighted_edgelist(padj_file_path):
        """
        makes a weighted edgelist from padj file, padj values are weights
        """

        try:
            edge_list = []
            with open(padj_file_path, "r") as padj_file:
                for line in padj_file:
                    line = line.split()
                    edge_list.append(line[0] + ":" + line[1] + "-" + line[2] + " " + line[3] + ":" + line[4] + "-" + line[5] + " " + line[6])

            return edge_list

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {padj_file}: {e}, PID: {os.getpid()}, TID: {tid}")
            raise

    @staticmethod
    def input_to_make_weighted_edgelist():
        """
        calls make_weighted_edgelist on all padj files
        """

        global first_print
        first_print = True

        padj_dir_path = SetDirectories.get_temp_dir() + "/padj"
        padj_dir = os.listdir(padj_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_temp_dir(), "weighted_edgelists")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        for file in padj_dir:
            full_path = os.path.join(padj_dir_path, file)
            if not os.path.isfile(full_path):
                print(f"File {file} does not exist.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            full_paths = [os.path.join(padj_dir_path, file) for file in padj_dir]
            futures = list(executor.map(Pipeline.make_weighted_edgelist, full_paths))
            for padj_file_path, future in zip(padj_dir, futures):
                try:
                    weighted_edgelist = future
                    output_filename = f"{os.path.basename(padj_file_path)[:-len('_padj.txt')]}_edgelist.txt"
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        for line in weighted_edgelist:
                            f.write(line + "\n")
                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {padj_file_path}: {e}, PID: {os.getpid()}, TID: {tid}")

    @staticmethod
    def make_edgelist(padj_file_path):

        """Makes edge list from padj file"""

        try:
            edge_list = []
            with open(padj_file_path, "r") as padj_file:
                for line in padj_file:
                    line = line.split()
                    edge_list.append(line[0] + ":" + line[1] + "-" + line[2] + "  " + line[3] + ":" + line[4] + "-" + line[5])

            custom_print(f"Finished processing file: {os.path.basename(padj_file_path)}, PID: {os.getpid()}, TID: {threading.get_ident()}")
            return edge_list

        except Exception as e:
            tid = threading.get_ident()
            print(f"Error in {padj_file}: {e}, PID: {os.getpid()}, tid: {tid}")
            raise

    @staticmethod
    def input_to_make_edgelist():
        """
        Calls make_edgelist on all padj files
        """

        global first_print
        first_print = True

        padj_dir_path = SetDirectories.get_temp_dir() + "/padj"
        padj_dir = os.listdir(padj_dir_path)

        # Create the output directory if it doesn't exist
        output_dir = os.path.join(SetDirectories.get_output_dir(), "edgelists")
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        else:
            shutil.rmtree(output_dir)
            os.mkdir(output_dir)

        full_paths = [os.path.join(padj_dir_path, file) for file in padj_dir]

        with concurrent.futures.ThreadPoolExecutor(max_workers=SetDirectories.get_threads()) as executor:
            futures = list(executor.map(Pipeline.make_edgelist, full_paths))
            for padj_file_path, future in zip(full_paths, futures):
                try:
                    edgelist = future
                    output_filename = f"{os.path.basename(padj_file_path)[:-len('_padj.txt')]}_edgelist.txt"
                    output_filepath = os.path.join(output_dir, output_filename)
                    with open(output_filepath, "w") as f:
                        for line in edgelist:
                            f.write(line + "\n")
                except Exception as e:
                    tid = threading.get_ident()
                    print(f"Error processing {padj_file_path}: {e}, PID: {os.getpid()}, TID: {tid}")