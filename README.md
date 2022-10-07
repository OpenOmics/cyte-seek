<div align="center">
   
  <h1>cyte-seek 🔬</h1>
  
  **_An Awesome Single-cell CITE-sequencing Pipeline_**

  [![Docker Pulls](https://img.shields.io/docker/pulls/skchronicles/chicyte)](https://hub.docker.com/repository/docker/skchronicles/chicyte) [![tests](https://github.com/OpenOmics/cyte-seek/workflows/tests/badge.svg)](https://github.com/OpenOmics/cyte-seek/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/cyte-seek/workflows/docs/badge.svg)](https://github.com/OpenOmics/cyte-seek/actions/workflows/docs.yml)<br>[![GitHub issues](https://img.shields.io/github/issues/OpenOmics/cyte-seek?color=brightgreen)](https://github.com/OpenOmics/cyte-seek/issues) [![GitHub license](https://img.shields.io/github/license/OpenOmics/cyte-seek)](https://github.com/OpenOmics/cyte-seek/blob/main/LICENSE)  
  
  <i>
    This is the home of the pipeline, cyte-seek. Its long-term goals: to accurately perform cell filtering, normalization, clustering, differential expression analysis, and cell type prediction like no pipeline before!
  </i>
</div>


## Overview
Welcome to cyte-seek! Before getting started, we highly recommend reading through [cyte-seek's documentation](https://openomics.github.io/cyte-seek/).

The **`./cyte-seek`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>cyte-seek <b>run</b></code>](https://openomics.github.io/cyte-seek/usage/run/): Run the cyte-seek pipeline with your input files.
 * [<code>cyte-seek <b>unlock</b></code>](https://openomics.github.io/cyte-seek/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>cyte-seek <b>cache</b></code>](https://openomics.github.io/cyte-seek/usage/cache/): Cache software containers locally.


**cyte-seek** is a comprehensive single-cell pipeline optimized for CITE-sequencing data. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from 10x sequencing technologies. The CITE-seq analysis pipeline starts from sample FASTQ files and performs initial processing of the data. Starting from Cell Ranger analysis, each sample undergoes cell filtering, data normalization, clustering, differential expression analysis, and cell type prediction. The processed samples are also integrated together to perform clustering, differential expression, and cell type prediction on a project level. Additionally, if genotype information is provided then genetic multiplexing of each sample is performed. The pipeline can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. 

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/cyte-seek/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/cyte-seek/issues). If you have any questions, please feel free to [start a discussion](https://github.com/OpenOmics/cyte-seek/discussions).

## Dependencies
**Requires:** `singularity>=3.5`  `snakemake>=6.0`

At the current moment, the pipeline uses a mixture of enviroment modules and docker images; however, this will be changing soon! In the very near future, the pipeline will only use docker images. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline will rely on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation

### Biowulf
Please clone this repository using the following commands:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/cyte-seek.git
cd cyte-seek/
# Add dependencies to $PATH
# Biowulf users should run
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module load snakemake singularity
# Get usage information
./cyte-seek -h
```

### LOCUS
Please clone this repository using the following commands:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/cyte-seek.git
cd cyte-seek/
# Add dependencies to $PATH
# LOCUS users should run
qrsh -l h_vmem=4G -pe threaded 4
module load snakemake/6.0.5-Python-3.9.2
# Download resource bundle (12G)
# in current working directory
wget https://hpc.nih.gov/~OpenOmics/cyte-seek/cyte-seek_v1.0.0.tar.gz
export cyte_refs="$PWD/references"
mkdir -p "$cyte_refs"
tar -xvf cyte-seek_v1.0.0.tar.gz -C "$cyte_refs"
# Update config to use locally
# downloaded reference files
sed -i "s@/data/OpenOmics/references@$cyte_refs@g" config/*.json
# Get usage information
./cyte-seek -h
```

## Contribute 

This site is a living document, created for and by members like you. cyte-seek is maintained by the members of OpenOmics and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/OpenOmics/cyte-seek).

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
