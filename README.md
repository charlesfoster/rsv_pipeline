# rsv_pipeline
A pipeline to facilitate genomic surveillance of RSV using the ARTIC amplicon scheme(s) and paired end short-read sequencing.

Author: Dr Charles Foster

# Starting out
## Clone the repo
To begin with, clone the 'main' branch of this github repository:

```
git clone https://github.com/charlesfoster/rsv_pipeline.git

cd rsv_pipeline
```

## Set up a conda environment
Ensure conda/mamba is installed. Easiest to use the Miniforge installer. See: https://conda-forge.org/download/

Next, install dependencies using `mamba`:

```
mamba env create -f environment.yml
```

# Running the pipeline
## Quick start
To run the pipeline, simply ensure you have changed to the pipeline directory, make sure the `config.yaml` file is accurate (see below), activate the conda environment, then run the Snakemake command.

```
cd rsv_pipeline
mamba activate rsv_pipeline
snakemake --cores 8
```

Choose an appropriate value for `--cores` based on your machine.

## Configuration
### Edit the config file
Within the main directory is a file called `config.yaml`. You need to make sure that its values correspond to what you want to analyse.

The file is divided into three broad sections: 

#### Input/Output:
* reads_dir: the directory where your reads for analysis are stored.
* outdir: the directory where your results should go.
* suffix1: the 'suffix' of your forward reads (see below)
* suffix2: the 'suffix' of your reverse reads (see below)
* remove_sample_index: how the sample should be parsed from filenames
* keep_trimmed_reads: should quality trimmed reads be kept?

The suffix parameters control how the pipeline identifies individual samples from a directory of reads. The easiest way to explain is with an example. By default, the pipeline assumes that reads are named exactly as they are after coming off a MiSeq/iSeq, i.e. `*_L001_R1_001.fastq.gz` (forward) and `*_L001_R1_001.fastq.gz` (reverse). Consider the following:

reads1 = `my_sample_S1_L001_R1_001.fastq.gz`

reads2 = `my_sample_S1_L001_R1_001.fastq.gz`

parsed sample name: `my_sample` if `remove_sample_index` == True, else `my_sample_S1` if `remove_sample_index` == False.

If you do not set the `suffix` parameters appropriately then no samples will be identified.

#### Analysis Parameters:
Various parameters for SNP calling and coverage thresholds. You will only occasionally need to edit these parameters.

#### Reference and Database Files:
Paths to reference and database files. You will normally not need to change these.



