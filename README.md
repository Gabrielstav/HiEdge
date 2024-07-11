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

#1. Clone the repository

```bash
git clone https://github.com/Gabrielstav/HiEdge.git
cd HiEdge
```

#2. Set up a virtal env

```bash
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
 ```

#3. Install dependancies

```bash
pip install -r requirements.txt
```

### Running the pipeline 

To run the pipeline, simply run the main.py script with the path to the configuration file as an argument (--config, -c flag).
If needed you can also name the run with the (--run_name, -r flag), if no name is given the run directory will be named after the current timestamp. 

```bash
python main.py -c path/to/config_file.yaml -r run_name
```



