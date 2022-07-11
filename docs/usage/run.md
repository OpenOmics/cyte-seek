# <code>cyte-seek <b>run</b></code>

## 1. About 
The `cyte-seek` executable is composed of several inter-related sub commands. Please see `cyte-seek -h` for all available options.

This part of the documentation describes options and concepts for <code>cyte-seek <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running cyte-seek pipeline. 

Setting up the cyte-seek pipeline is fast and easy! In its most basic form, <code>cyte-seek <b>run</b></code> only has *five required inputs*.

## 2. Synopsis
```text
$ cyte-seek run [--help] [--mode {slurm,uge,local}] [--job-name JOB_NAME] \
      [ --dry-run] [--silent] [--sif-cache SIF_CACHE] \
      [--singularity-cache SINGULARITY_CACHE]  \
      [--tmpdir TMP_DIR] [--threads THREADS] \
      [--pre-rna] [--force-cells] \
      [--num-cells NUM_CELLS] \
       --input INPUT [INPUT ...] \
       --output OUTPUT \
       --genome {hg38, ...} \
       --libraries LIBRARIES \
       --features FEATURES
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument, an output directory to store results via `--output` argument, a reference genome for alignment and cell-type prediction via the `--genome` arguement, a [libraries file](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis) via the `--libraries` argument, and a [features barcode](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis) file via a `--features` argument. 

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ file(s).**  
> *type: file(s)*  
> 
> One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should seperated by a space. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
> 
> ***Example:*** `--input .tests/*.R?.fastq.gz`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/output`

---  
  `--genome {hg38, ...}`
> **Reference genome.**   
> *type: string*
>   
> This option defines the reference genome of the samples. cyte-seek does comes bundled with prebuilt reference files for human and mouse samples, e.g. hg38. Please select one of the following options: hg38. Please note that the mouse reference genome, mm10, is coming soon! 
> 
> ***Example:*** `--genome hg38`

---  
  `--libraries LIBRARIES`
> **Libraries file.**   
> *type: file*
>   
> A CSV file containing information about each library. It contains each sample's name, flowcell, demultiplexed name, and library type. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).

> *Here is an example libraries.csv file:*
> ```
> Name,Flowcell,Sample,Type
> IL15_LNs,H7CNNBGXG,IL15_LNs,Gene Expression
> IL15_LNs,H7CT7BGXG,IL15_LNs_BC,Antibody Capture
> ``` 

> *Where:* 

> - *Name:* name of the sample passed to CellRanger.  
> - *Flowcell:* The flowcell ID that contains the FASTQ files for this set of data.  
> - *Sample:* Name that was used when demultiplexing, this should match the FASTQ files.  
> - *Type:* library type for each sample. List of supported options:  
>        * Gene Expression
>        * CRISPR Guide Capture
>        * Antibody Capture
>        * Custom

> ***Example:*** `--libraries libraries.csv`


---  
  `--features FEATURES`
> **Features file.**   
> *type: file*
>   
> A feature reference CSV file containing information for processing a feature barcode data. This file should contain a unique ID for the feature, a human readable name, sequence, feature type, read, and pattern. More information about the libraries file and its requirements can be found on the [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).

> *Here is an example features.csv file:*
> ```
> id,name,sequence,feature_type,read,pattern
> CITE_CD64,CD64,AGTGGG,Antibody Capture,R2,5PNN(BC)
> CITE_CD8,CD8,TCACCGT,Antibody Capture,R2,5PNNN(BC)
> ``` 

> *Where:*  

> - *id:* Unique ID for this feature. Must not contain whitespace, quote or comma characters. Each ID must be unique and must not collide with a gene identifier from the transcriptome. 
> - *name:* Human-readable name for this feature. Must not contain whitespace. 
> - *sequence:* Nucleotide barcode sequence associated with this feature, e.g. the antibody barcode or sgRNA protospacer sequence. 
> - *read:* Specifies which RNA sequencing read contains the Feature Barcode sequence. Must be R1 or R2, but in most cases R2 is the correct read. 
> - *pattern:* Specifies how to extract the sequence of the feature barcode from the read.

> - *Type:* Type of the feature. List of supported options:  
>        * Gene Expression
>        * CRISPR Guide Capture
>        * Antibody Capture
>        * Custom

> ***Example:*** `--features features.csv`

### 2.2 Analysis options

Each of the following arguments are optional, and do not need to be provided. 

  `--num-cells NUM_CELLS`            
> **Expected number of recovered cells.**  
> *type: int*  
> *default: 3000*
>
> Overrides the expected number of cells passed to cellranger.
>
> ***Example:*** `--num-cells 4000`

---  
  `--force-cells`            
> **Use expected cells counts.**  
> *type: boolean*
> 
> Force pipeline to use the expected number of cells, which will bypass the cell detection algorithm.
>
> ***Example:*** `--force-cells`


---  
  `--pre-mrna`            
> **Consider pre-mRNA.**  
> *type: boolean*
> 
> Retain intronic sequences for consideration of pre-mRNA.
>
> ***Example:*** `--pre-mrna`

---  
  `--demuxlet`            
> **Perform demuxlet analysis.**  
> *type: boolean*
> 
> Perform demuxlet analysis for this project. This option requires a vcf file is provided. Please see the option below for more information.
>
> ***Example:*** `--demuxlet`

---  
  `--vcf VCF`            
> **VCF file used for demuxlet analysis.**  
> *type: file*  
> 
> This option should be used with the demuxlet option above.
>
> ***Example:*** `--vcf analysis.vcf`


---  
  `--patient-list PATIENT_LIST`            
> **Define patients associated with each sample for demuxlet analysis.**  
> *type: file*
> 
> A CSV file used to define the patients associated with each single cell sample used in demuxlet analysis.

> *Here is an example features.csv file:*
> ```
> Sample1,Sample_2
> patient1,patient2
> patient2,patient4
> patient3,
> ``` 
>
> *Where:*  

> - *first row*: contains the single cell sample names matching the Cell Ranger output
> - *following row(s):* contain the patients associated with each single cell sample, using the same patient ID listing as listed in the vcf file.

> ***Example:*** `--patient-list demuxlet.csv`

### 2.4 Orchestration options

Each of the following arguments are optional, and do not need to be provided. 

  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
> 
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,uge,local}`  
> **Execution Method.**  
> *type: string*  
> *default: slurm*
> 
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm, uge, or local. It is recommended to run cyte-seek on a cluster to reduce run times. At the current moment, the pipeline supports the following job schedulers: SLURM, UGE.
> 
> ***slurm***    
> The slurm execution backend will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). This method will submit jobs to a cluster using sbatch. Please set the mode to slurm when running the pipeline on Biowulf. 
>
> ***uge***    
> The uge execution backend will submit jobs to the [UGE workload manager](https://en.wikipedia.org/wiki/Univa_Grid_Engine). This method will submit jobs to a cluster using qsub. Please set the mode to uge when running the pipeline on LOCUS. 
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a local execution mode. We do not recommend using this option.
> 
> ***Example:*** `--mode uge`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*
> *default: cyte-seek*
> 
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "cyte-seek".
> 
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `--output OUTPUT/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
> 
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `cyte-seek cache` subcommand can be used to create a local SIF cache. Please see `cyte-seek cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running cyte-seek with this option when ever possible.
> 
> ***Example:*** `--singularity-cache /data/$USER/SIFs`

---  
  `--threads THREADS`   
> **Max number of threads for local processes.**  
> *type: int*  
> *default: 2*
> 
> Max number of threads for local processes. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
> 
> ***Example:*** `--threads 12`


---  
  `--tmp-dir TMP_DIR`   
> **Max number of threads for each process.**  
> *type: path*  
> *default: `/tmp/$USER`*
> 
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOBID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
> 
> ***Example:*** `--tmp-dir /hpcdata/scratch/$USER`

### 2.5 Miscellaneous options  
Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example

### 3.1 Biowulf
```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./cyte-seek run --input .tests/*.R?.fastq.gz \
    --output /data/$USER/cyte-seek_output \
    --features .tests/features.csv \
    --libraries libraries.csv \
    --genome hg38 \
    --mode slurm \
    --tmp-dir '/lscratch/$SLURM_JOBID' \
    --dry-run


# Step 2B.) Run the cyte-seek pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode on Biowulf.
./cyte-seek run --input .tests/*.R?.fastq.gz \
    --output /data/$USER/cyte-seek_output \
    --features .tests/features.csv \
    --libraries libraries.csv \
    --genome hg38 \
    --mode slurm \
    --tmp-dir '/lscratch/$SLURM_JOBID'
```

### 3.2 LOCUS
```bash 
# Step 1.) Grab an interactive node,
# do not run on head node!
qrsh -l h_vmem=4G -pe threaded 4
module load snakemake/6.0.5-Python-3.9.2

# Step 2A.) Dry-run the pipeline
/hpcdata/rtb/NCBR/cyte-seek/v1.0.0/cyte-seek run \
    --input .tests/*.R?.fastq.gz \
    --output /hpcdata/scratch/$USER/cyte-seek_output \
    --features .tests/features.csv \
    --libraries libraries.csv \
    --genome hg38 \
    --mode uge \
    --dry-run

# Step 2B.) Run the cyte-seek pipeline
# The slurm mode will submit jobs to 
# the cluster. It is recommended running 
# the pipeline in this mode on Biowulf.
/hpcdata/rtb/NCBR/cyte-seek/v1.0.0/cyte-seek run \
    --input .tests/*.R?.fastq.gz \
    --output /hpcdata/scratch/$USER/cyte-seek_output \
    --features .tests/features.csv \
    --libraries libraries.csv \
    --genome hg38 \
    --mode uge 
```