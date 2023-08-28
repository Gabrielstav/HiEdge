# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.

import argparse
import pathlib as pl

def arg_parser():

    help_message = """
    Hi-C data pipeline for processing Hi-C data from HiC-Pro to statistically significant edge-lists.

    Essential arguments include the path to the YAML configuration file and paths for the input and output directories.
    Please refer to the comprehensive documentation for details about the configuration settings in the YAML file.
    https://github.com/Gabrielstav/HiEdge

    Example Usage:
    python main_script.py -c /path/to/config.yaml -i /path/to/input -o /path/to/output
    """

    parser = argparse.ArgumentParser(description=help_message, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-c", "--config", help="Path to the YAML configuration file.", required=True)
    parser.add_argument("-i", "--input_dir", help="Path to the input directory containing HiC-Pro output folders (bed and matrix files).", required=True)
    parser.add_argument("-o", "--output_dir", help="Path to the output directory to store processed data.", required=True)

    args = parser.parse_args()

    # Ensure config file exists
    if not pl.Path(args.config).exists():
        parser.error("The provided YAML configuration file doesn't exist!")

    return args

help_message = "HiC_Pipeline for processing Hi-C data from HiC-Pro to statistically significant edge-lists. \n\n" \
               "INPUT DIR: -i\n" \
               "Directory containing HiC-Pro output folders (bed and matrix files) is set as input directory. Any folder can be the input, as long as it contains the HiC-Pro output folders (raw, matrix) for one HiC-Pro run. \n\n" \
               "OUTPUT DIR: -o \n" \
               "Any directory to output processed data is set as output directory.  \n\n" \
               "REFERENCE DIR: -r \n" \
               "Directory containing reference genome files is set as reference directory. Put the reference files in a folder called hg19. \n" \
               "This directory should contain the following files: \n" \
               "    cytoBand_hg19.txt (USCS cytoband reference file: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/) \n" \
               "    hg19-blacklist.v2.bed (Encode hg19 blacklisted regions: https://github.com/Boyle-Lab/Blacklist/tree/master/lists) \n\n" \
               "NCHG PATH: -n \n" \
               "Path to NCHG executable is set as NCHG path. NCHG is a statistical tool written in C++ for calculating p-values for Hi-C interactions using the NCHG distribution. \n" \
               "It can be found here: https://github.com/Chrom3D/preprocess_scripts/blob/master/NCHG_hic.zip.\n" \
               "Credit: Jonas Paulsen (2017). \n\n" \
               "NORM OPTION: -m \n" \
               "Specifies if normalized data or raw data is processed to edge lists. Options: raw, iced, norm, normalized. \n" \
               "Raw data is default, since the NCHG tool uses raw counts to calculate significance. \n" \
               "If raw is selected, the script will look for raw data in the HiC-Pro output folder. If iced is selected, the script will look for ICE normalized data in the HiC-Pro output folder. \n" \
               "\n\n" \
               "If no arguments are given, the script will run with the hardcoded paths set in the SetDirectories class. Meaning, it's possible to run without providing arguments. \n" \
               "For instance, set the NCHG and reference paths as hardcoded, and provide input and output directories for each run. \n\n" \
               "INTER-CHROMOSOMAL INTERACTIONS: -inter. \n" \
               "If specified, run NCHG with the -i flag, meaning interchromosmal interactions will be included in calculating significance. \n\n" \
               "INTRA-CHROMOSOMAL INTERACTIONS: -intra. \n" \
               "If specified, run NCHG without the -i flag, meaning only intrachromosmal interactions will be used to calcualte significance. This is the default. \n\n" \
               "MIXED INTERACTIONS: -mixed \n" \
               "If this flag is set, the NCHG script will consider if inter- or intrachromosomal interactions are appropriate per input file. \n" \
               "Allows for running mixed input files, some of which contain inter-chromosomal interactions and some that do not. \n\n" \
               "THREADS: -t \n" \
               "Int: Number of threads to use for processing. Default is cores available on machine. \n\n" \
               "FDR THRESHOLD: -f \n" \
               "Float: FDR threshold for significance. Default is 0.05. \n\n" \
               "RESOLUTIONS: -res \n" \
               "INT: Resolution values can be provided to run the pipeline on specific resolutions (needs to match resolutions available in HiC-Pro output folder). \n" \
               "The script will run on all resolutions if no resolution values are provided. \n\n" \
               "EXECUTOR: -e \n" \
               "String: Executor to use for parallelization. Options: thread, process. Default is multiprocessing. \n\n" \
               "NO_SPLIT: -ns \n" \
               "If specified, do not split input files to NCHG by each chromosome. Default is to split files, to allow increased multiprocessing of NCHG. \n\n" \
               "N_QUANTILES: -nq \n" \
               "Int: Number of quantiles to use in NCHG for estimation of expected number of interactions given a genomic distance. Default is 100. \n\n" \

    parser = argparse.ArgumentParser(description=help_message, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input_dir", help="Directory containing HiC-Pro output folders (bed and matrix files)", required=False)
    parser.add_argument("-o", "--output_dir", help="Any directory to output processed data", required=False)
    parser.add_argument("-r", "--reference_dir", help="Directory containing reference genome files.", required=False)
    parser.add_argument("-n", "--nchg_path", help="Path to NCHG executable", required=False)
    parser.add_argument("-m", "--norm_option", help="Normalization option, default is to use raw counts.", choices=["raw", "iced", "norm", "normalized"], required=False)
    parser.add_argument("-inter", "--interchromosomal_interactions", help="Consider inter-chromosomal interactions to find significant interactions.", action="store_true", required=False, default=None)
    parser.add_argument("-intra", "--intrachromosomal_interactions", help="Consider only intra-chromosomal interactions to find significant interactions (default).", action="store_true", required=False, default=None)
    parser.add_argument("-mixed", "--mixed_interactions", help="Consider inter-chromosomal and intra-chromosmal interactions on a file-by-file basis.", action="store_true", required=False, default=None)
    parser.add_argument("-t", "--threads", help="Int: Number of threads to use for processing. Default is cores available on machine. Always specify on HPC cluster.", required=False)
    parser.add_argument("-f", "--fdr_threshold", help="Float: FDR threshold for significance. Default is 0.05.", required=False, type=float)
    parser.add_argument("-res", "--resolutions", help="Int: Resolution values can be provided to run the pipeline on specific resolutions.", required=False, nargs="+", type=int)
    parser.add_argument("-e", "--executor", choices=["m", "t", "multi", "mp", "th", "thread", "processing", "multiprocessing", "threading"], default="multiprocessing", help="Choose between multiprocessing and threading for the NCHG script method. Default is multiprocessing.")
    parser.add_argument("-ns", "--no_split", action="store_true", required=False, help="Do not split input files to NCHG by chromosomes.")
    parser.add_argument("-nq", "--n_quantiles", help="Int: Number of quantiles to use for NCHG. Default is 100.", required=False, type=int)
    args = parser.parse_args()

    # Sets args to None if not provided in command line
    input_directory = args.input_dir
    output_directory = args.output_dir
    reference_directory = args.reference_dir
    nchg_executable_path = args.nchg_path
    norm_option = args.norm_option
    interchromosomal_interactions = args.interchromosomal_interactions
    intrachromosomal_interactions = args.intrachromosomal_interactions
    mixed_interactions = args.mixed_interactions
    fdr_threshold = args.fdr_threshold
    threads = args.threads
    resolutions = args.resolutions
    executor_type = args.executor
    no_split = args.no_split
    n_quantiles = args.n_quantiles

    args = parser.parse_args()

    return args