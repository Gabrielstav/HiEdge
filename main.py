# Copyright Gabriel B. Stav. Licensed under the terms of the Apache 2.0 license. See LICENSE in the project root.


# Import modules
import concurrent.futures
import shutil
import pybedtools as pbt
import os as os
import subprocess as sp
import math
from statsmodels.sandbox.stats import multicomp
import time as time
from tqdm import tqdm
import re as re
import argparse as argparse
import pickle as pickle
import threading as threading
from collections import defaultdict


# Call each run method from each class in the filtering,
# first filtering, then processing, then statistical testing then plotting/output


def run_pipeline():
    """
    Call selected methods of the HiC_Pipeline, in the order specified
    """

    global first_print
    first_print = True

    start_time = time.time()

    # List of methods to call
    method_names = [
        (lambda: Pipeline.input_to_make_bedpe(Pipeline_Input.group_files(SetDirectories.get_input_dir()))),
        "input_to_remove_blacklist",
        "input_to_remove_cytobands",
        "input_to_nchg",
        "input_to_adjust_pvalues",
        "input_to_make_edgelist"
    ]

    # Call each method once
    for method in tqdm(method_names):
        if type(method) == str:  # Check if method name is a string
            method = getattr(Pipeline, method)  # Get the method reference
        method()  # Call the method

    # Print runtime on completion
    end_time = time.time()
    print(f"HiC_Pipeline completed in {end_time - start_time:.2f} seconds. ({(end_time - start_time) / 60:.2f} minutes, {((end_time - start_time) / 60) / 60:.2f} hours).")


if __name__ == "__main__":
    run_pipeline()




