<div align="center">
   
  <h1 style="font-size: 250%">cyte-seek 🔬</h1>  
  
  <b><i>An Awesome Single-cell CITE-sequencing Pipeline</i></b><br>
  <a href="https://hub.docker.com/repository/docker/skchronicles/chicyte">
    <img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/skchronicles/chicyte">
  </a>  
  <a href="https://github.com/OpenOmics/cyte-seek/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/cyte-seek/workflows/tests/badge.svg">
  </a>  
  <a href="https://github.com/OpenOmics/cyte-seek/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/cyte-seek/workflows/docs/badge.svg">
  </a><br>  
  <a href="https://github.com/OpenOmics/cyte-seek/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/cyte-seek?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/cyte-seek/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/cyte-seek">
  </a>
  
  <p>
    This is the home of the pipeline, cyte-seek. Its long-term goals: to accurately perform cell filtering, normalization, clustering, differential expression analysis, and cell type prediction like no pipeline before!
  </p>

</div> 


## Overview
Welcome to cyte-seek's documentation! This guide is the main source of documentation for users that are getting started with the [single-cell CITE-sequencing pipeline](https://github.com/OpenOmics/cyte-seek/). Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

The **`./cyte-seek`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>cyte-seek <b>run</b></code>](usage/run.md): Run the cyte-seek pipeline with your input files.
 * [<code>cyte-seek <b>unlock</b></code>](usage/unlock.md): Unlocks a previous runs output directory.
 * [<code>cyte-seek <b>cache</b></code>](usage/cache.md): Cache remote resources locally, coming soon!

**cyte-seek** is a comprehensive single-cell pipeline optimized for CITE-sequencing data. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from 10x sequencing technologies. The CITE-seq analysis pipeline starts from sample FASTQ files and performs initial processing of the data. Starting from Cell Ranger analysis, each sample undergoes cell filtering, data normalization, clustering, differential expression analysis, and cell type prediction. The processed samples are also integrated together to perform clustering, differential expression, and cell type prediction on a project level. Additionally, if genotype information is provided then genetic multiplexing of each sample is performed. The pipeline can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. 

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/cyte-seek/issues). If you have any questions, please feel free to [start a discussion](https://github.com/OpenOmics/cyte-seek/discussions).

## Contribute 

This site is a living document, created for and by members like you. cyte-seek is maintained by the members of OpenOmics and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/cyte-seek).

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
