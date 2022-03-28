# Functions and rules for processing CITE-seq data

# Function defitions
def count_premrna(wildcards):
    """
    Wrapper to decide whether to include introns for counting.
    See config['options']['pre_mrna'] for the encoded value.   
    """
    if pre_mrna:
        return('--include-introns')
    else:
        return('')

def count_expect_force(wildcards):
    """
    Wrapper to decide whether to force the expected number of cells.
    See config['options']['force_cells'] for the encoded value.
    """
    if force_cells:
        return('--force-cells')
    else:
        return('--expect-cells')


# Rule definitions
rule librariesCSV:
    output: 
        expand(join(workpath, "{sample}_libraries.csv"), sample=samples)
    params: 
        rname = "libcsv",
        fastq = ",".join(input_paths_set),
        libraries = libraries,
        create_libs = join("workflow", "scripts", "create_library_files.py"),
    shell: 
        """
        python {params.create_libs} \\
            {params.libraries} \\
            {params.fastq}
        """


rule count:
    input: 
        lib = join(workpath, "{sample}_libraries.csv"),
        features = features
    output: 
        join(workpath, "{sample}", "outs", "web_summary.html")
    log: 
        err = "run_{sample}_10x_cellranger_count.err", 
        log ="run_{sample}_10x_cellranger_count.log"
    params: 
        rname = "count",
        batch = "-l nodes=1:ppn=16,mem=96gb", 
        prefix = "{sample}", 
        numcells = lambda wildcards:s2c[wildcards.sample], 
        transcriptome = config["references"][genome]["transcriptome"],
        premrna = count_premrna, 
        cells_flag = count_expect_force
    envmodules:
        config["tools"]["cellranger"]
    shell: 
        """
        # Remove output directory
        # prior to running cellranger
        if [ -d '{params.prefix}' ]; then
            rm -rf '{params.prefix}/'
        fi

        cellranger count \\
            --id={params.prefix} \\
            {params.cells_flag}={params.numcells} \\
            --transcriptome={params.transcriptome} \\
            --libraries={input.lib} \\
            --feature-ref={input.features} \\
            {params.premrna} \\
        2>{log.err} 1>{log.log}
        """


rule summaryFiles:
    input: 
        expand(join(workpath, "{sample}", "outs", "web_summary.html"), sample=samples)
    output: 
        join(workpath, "finalreport", "metric_summary.xlsx"),
        expand(join(workpath, "finalreport", "summaries", "{sample}_web_summary.html"), sample=samples)
    params: 
        rname = "sumfile",
        batch = "-l nodes=1:ppn=1",
        summarize = join("workflow", "scripts", "generateSummaryFiles.py"),
    shell: 
        """
        python2 {params.summarize}
        """


rule aggregateCSV:
    input: 
        expand(join(workpath, "{sample}", "outs", "web_summary.html"), sample=samples)
    output: 
        join(workpath, "AggregatedDatasets.csv"),
    params: 
        rname = "aggcsv",
        batch = "-l nodes=1:ppn=1",
        outdir = workpath,
        aggregate = join("workflow", "scripts", "generateAggregateCSV.py"),
    shell: 
        """
        python2 {params.aggregate} {params.outdir}
        """


rule aggregate:
    input: 
        csv=join(workpath, "AggregatedDatasets.csv"),
    output: 
        touch(join(workpath, 'aggregate.complete')),
    log: 
        err="run_10x_aggregate.err", 
        log="run_10x_aggregate.log",
    params: 
        rname = "agg",
        batch = "-l nodes=1:ppn=16,mem=96gb",
    envmodules:
        config["tools"]["cellranger"]
    shell: 
        """
        cellranger aggr \\
            --id=AggregatedDatasets \\
            --csv={input.csv} \\
            --normalize=mapped \\
        2>{log.err} 1>{log.log}
        """
