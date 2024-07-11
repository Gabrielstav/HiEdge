# HiEdge

## Overview

HiEdge is a Hi-C pipeline used to easily and automatically process output data from Hi-C Pro to significant interactions. All you need to do it provide the input directory containing the Hi-C data, and the pipeline does the rest.

Supported input formats are the HiC-Pro output files (.BED, .matrix, .matrix.balanced & .bias). 
The pipeline matches files from the same resolution and experiment and aggregates them to one datastructure.
This allows for easy processing of multiple resolutions and experiments in parallel - and per-dataset filtering options - including processing of inter- and intrachromosomal interactions in the same run. 
These filtering options include chromosomes, specific genomic regions, blacklisted regions and ceontromeres. 
This means that you can easily filter out unwanted data and only keep the significant interactions on a per-dataset and per-resolution basis. 
Each filtering option can be customized or turned entirely off - allowing for modular and flexible processing of the data.

After filtering the data, singificance testing is done - which is handled separately for inter- and intrachromosomal interactions due to the differences in background noise and expected interactions.
For intrachromosomal data the bins are aggregated to "metabins" and these metabins are then used to generate a monotonically decreasing spline fit to the genomic distance over interaction frequency.
This spline is then used to generate a null model for the data, which is used to calculate the p-values using a binomial survival test, before doing FDR correction.
For interchromosomal data there is no distance-dependant decay, instead the interaction frequencies are used directly in the survival test before FDR correction.
The end result is a list of significant interactions, with p-values and adjusted p-values, that can be used for further analysis.
The output format can be specified as needed. 

All settings are specified in a configuration file, which allows for easy customization of the pipeline.
One run of the pipeline can contain multiple datasets, and each dataset can contain multiple resolutions and experiments - but one run of the pipeline corresponds to one configuration file.


![1720685919072](HiEdge_fig.png)

### Installation

1. Clone the repository

```bash
git clone https://github.com/Gabrielstav/HiEdge.git
cd HiEdge
```

2. Set up a virtual environment

```bash
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
 ```

3. Install dependencies

```bash
pip install -r requirements.txt
```

### Running the pipeline 

To run the pipeline, simply run the main.py script with the path to the configuration file as an argument (--config, -c flag).
If needed you can also name the run with the (--run_name, -r flag), if no name is given the run directory will be named after the current timestamp. 

```bash
python main.py -c path/to/config_file.yaml -r run_name
```

### Configuration 

The configuration file is a YAML file that specifies the settings for the pipeline.

# Paths
- input_dir: Directory containing HiC-Pro output files
- run_dir: Output directory containing HiEdge results
- output_dir: Output directory for HiEdge results
- temp_dir: Temporary directory for storing intermediate files

# Reference data
- blacklist_dir: Directory containing blacklisted regions (must be in ENCODE BED format and correspond to the genome version used in the HiC-Pro experiment)
- cytoband_dir: Directory containing cytoband regions (must be in UCSC BED format and correspond to the genome version used in the HiC-Pro experiment)

The current version only has support for hg19, but hg38 support is planned for the future.

# Pipeline settings
- reference_genome: Reference genome used in HiC-Pro experiment upstream of pipeline, used for filtering blacklisted and cytogenic regions (hg19 or hg38).
- hicpro_raw_dirname: Directory containing raw HiC-Pro matrix output (dafault name by HiC-Pro is always: raw).
- hicpro_norm_dirname: Directory containing balanced HiC-Pro matrix output and bias matrices (dafault name by HiC-Pro is always: iced).
- output_type: Output type for HiEdge results (default, verbose, qval). This denotes how many metrics are outputted in the final results. 
  Default contains the BED coordinates with the p-values, interaction counts and adjusted p-values, verbose contains all metrics and qval only contains the BED format with adjusted p-values.
- output_format: Output format for HiEdge results (csv, hdf5, parquet, txt).
- make_plots: If true, make plots for HiEdge results (spline fit, distance dependant decay, p-value distribution, q-value distribution, etc). 

# Execution settings
- interaction_type: Type of interactions to consider (inter, intra, mixed). Mixed will process all interactions, but will perform separate statistical testing for inter and intra interactions.
- iced_data: If set to true, use balanced matrix files by iterative correction (ic) from Hi-C Pro (iced matrices) instead of raw matrix files (not recommended, statistical model assumes raw matrices).
- round_iced_matrices: If true, round iced matrices to integer values.
- intra_resolutions: List of integers of resolutions in base pairs for intrachromosomal data to be processed (if empty, process all resolutions for which data is available).
- inter_resolutions: List of integers of resolutions in base pairs to consider for interchromosomal data (if empty, process all resolutions for which data is available).

# Filtering settings
- filter_blacklist: If true, filter out interactions in blacklisted regions (blacklist coordinate file).
- filter_cytobands: If true, filter out interactions in cytoband regions (specified by cytoband coordinate file).
- remove_chromosomes: List of strings of chromosome(s) to be removed from dataset. Empty list means no chromosomes are selected.
- select_chromosomes: List of strings of chromosome(s) to include in dataset, all chromosomes not selected are omitted. Empty list means no specific chromosomes are selected.
- filter_self_interactions: If true, remove self-interacting bins (interactions between hic-contacts in the same bin - cycles). 
- select_specific_regions: If true, select specific regions for filtering. Must be set to true if select_regions is specified.
- select_regions: Dictionary of chromosome(s) and list of regions to include in dataset. Regions are specified as strings in the format "start-end". 
- omit_specific_regions: If true, omit specific regions for filtering. Must be set to true if omit_regions is specified.
- omit_regions: Dictionary of chromosome(s) and list of regions to exclude from dataset. Regions are specified as strings in the format "start-end".
- use_interaction_distance_filters: If true, filter out interactions based on distance. Must be set to true if interaction_distance_filters is specified.
- interaction_distance_filters: Dictionary of distance thresholds and minimum and maximum distances for interactions to be considered. 
  Each key is a resolution, where min and max distance are the range withtin which interactions are considered.

# Statistical settings
- spline_passes: Number of spline passes for spline fitting - increasing the number of spline passes can improve the fit but also increase the runtime - default is 1 (additonal passes not yet supported).
- fdr_threshold: False Discovery Rate threshold for adjusting p-values by Benjamini-Hochberg procedure (default is 0.05). 
- metabin_occupancy: Number of metabins to use for occupancy calculation (number of bins used in spline fit - default is 200). 
- use_hicpro_bias: If true, use bias files from HiC-Pro for adjusting expected interaction frequencies before statistical testing. If false, no bias files are used. 
  Using bias files for normalization is recommended, as it can improve the fit of the spline model at the cost of slight increase in runtime.
- bias_lower_bound: Lower bound for bias values for filtering (default is 0.5).
- bias_upper_bound: Upper bound for bias values for filtering (default is 2).
- use_filtered_data_for_average_contact_probability: If true, use filtered data for average contact probability calculation used in spline fitting. If false, use raw (unfiltered) data.
- use_sequential_fdr: If true, use sequential FDR correction for multiple testing. If false, use partition-based FDR correction.

