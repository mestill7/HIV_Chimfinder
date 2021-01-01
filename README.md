# HIV_chim:

This repository hosts an pipeline for detecting transcriptional chimeras between a mammalian host genome (e.g. human or mouse) and a virus (e.g. HIV), created using Snakemake. Designed for paired-end RNA-seq. 

## Dependencies:
- [Anaconda](https://conda.io/docs/user-guide/install/linux.html) 

- BBMap (https://sourceforge.net/projects/bbmap/)

## Installation:
Clone this repository and change into the cloned NGS-Data-Charmer directory. 

To create an environment using the environment.yaml file, type the following:

`conda env create -f environment.yaml`

This will create a conda environment called "hiv_chim".

You will also have to download and unzip BBMap to a directory (for example, ~/Tools). 
You can test the bbmap install with the following command (modify paths as necessary):
~/Tools/bbmap/filterbyname.sh --help

## Usage note:

You must manually activate the conda environment prior to running the sh files. Type the following to activate the environment:

`conda activate hiv_chim

export PATH=$PATH:~/Tools/bbmap`

## Usage on an LSF cluster (highly recommended):

Copy the config.yaml, run\_snakemake\_cluster.sh, cluster.json and Snakefile to your NGS project directory. This directory should also contain a directory called 'fastq' wherein all fastq files are placed. Make sure the project directory structure is as follows:
```
.
├── cluster.json
├── config.yaml
├── fastq
│   ├── SampleA_R1_001.fastq.gz
│   └── SampleA_R2_001.fastq.gz
├── run_snakemake_cluster.sh
└── Snakefile (required for all analyses)
```
Make the required changes to the config.yaml and cluster.json file.

Next, type `nohup sh run_snakemake_cluster.sh &` (to run in background).

## Usage on a local machine:

Copy the config.yaml, run\_snakemake.sh and Snakefile to your NGS project directory. This directory should also contain a directory called 'fastq' wherein all the fastq files are placed. Make sure the project directory structure is as follows:
```
.
├── config.yaml
├── fastq
│   ├── SampleA_R1_001.fastq.gz
│   └── SampleA_R2_001.fastq.gz
├── run_snakemake.sh
└── Snakefile
```
Make the required changes to the config.yaml file.

Finally, type `sh run_snakemake.sh` followed by the maximum number of CPU cores to be used by snakemake. For example, type `sh run_snakemake.sh 2` for 2 CPU cores. You can also type `nohup sh run_snakemake.sh 2 &` to run the pipeline in background.

## Steps in pipeline:

![ScreenShot](/dag/dag.png)

## Examining results for transcriptional chimeras:
Your "output/spades_\*/transcripts.sam" will contain the alignments of the *de novo* transcriptome to the host genome (e.g. human) supplemented with a chromosome corresponding to the viral genome (e.g. HIV). To find chimeric reads, you can search for chimeras with the "XT:" flag.

`grep "XT:" output/spades_*/transcripts.sam | grep "HIV"`

## File naming requirements

It is best practice to ensure that sequencing files are properly paired (if paired-end sequencing was performed) and are named with consistent fastq file endings (e.g. ".R1.fq.gz" and ".R2.fq.gz") prior to analyzing sequencing reads. However, the Snakemake pipeline is robust to mixtures of fastq file endings. It does this by detecting the most common forward read file ending (e.g. "\*.R1.fq.gz"), then renaming any files that do not conform to the most common fastq file ending.

Please note that input fastq file names that fail to conform to any of the expected fastq filename endings (for example, "\_R1\_001.fastq.gz",".R1.fq.gz", "fq.gz", "fastq.gz", "fastq", and "fq" are all examples of allowed fastq file endings) will be ignored by Snakemake (e.g. files named "Treatment.gz" or "Treatment.txt" or "Treatment" will be ignored). Mixtures of ".gz" and non-gzipped *fastq* files are NOT allowed, and will result in an error. If you do have such a mixture, please gzip the uncompressed fastq files with the "gzip" command OR unzip any compressed fastq files prior to running the pipeline. 
