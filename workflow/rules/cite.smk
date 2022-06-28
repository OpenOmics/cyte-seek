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

def demuxlet_analysis(wildcards):
    if demuxlet:
        return(join(workpath, 'demuxlet', 'output', f'{wildcards.sample}', f'{wildcards.sample}.best'))
    else:
        return()

def seurat_optional_params(wildcards):
    if demuxlet:
        return(join(workpath, 'demuxlet', 'output', f'{wildcards.sample}', f'{wildcards.sample}.best'))
    else:
        return('')

def demuxlet_input(wildcards):
    if patient_list != 'None':
        return(join(workpath, 'demuxlet', 'output', f'{wildcards.sample}', 'patient_list'))
    else:
        return()

def demuxlet_flag(wildcards):
    if patient_list != 'None':
        return( "--sm-list " + join(workpath, 'demuxlet', 'output', f'{wildcards.sample}', 'patient_list'))
    else:
        return('')


# Rule definitions
rule librariesCSV:
    output:
        expand(join(workpath, "{sample}_libraries.csv"), sample=lib_samples)
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
        expand(join(workpath, "{sample}", "outs", "web_summary.html"), sample=lib_samples)
    output:
        join(workpath, "finalreport", "metric_summary.xlsx"),
        expand(join(workpath, "finalreport", "summaries", "{sample}_web_summary.html"), sample=lib_samples)
    params:
        rname = "sumfile",
        batch = "-l nodes=1:ppn=1",
        summarize = join("workflow", "scripts", "generateSummaryFiles.py"),
    shell:
        """
        python {params.summarize}
        """


rule aggregateCSV:
    input:
        expand(join(workpath, "{sample}", "outs", "web_summary.html"), sample=lib_samples)
    output:
        join(workpath, "AggregatedDatasets.csv"),
    params:
        rname = "aggcsv",
        batch = "-l nodes=1:ppn=1",
        outdir = workpath,
        aggregate = join("workflow", "scripts", "generateAggregateCSV.py"),
    shell:
        """
        python {params.aggregate} {params.outdir}
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

rule seurat:
    input:
        join(workpath, "{sample}", "outs", "web_summary.html"),
        demuxlet_analysis
    output:
        rds = join(workpath, "seurat", "{sample}", "seur_cite_cluster.rds")
    log:
        join("seurat", "{sample}", "seurat.log")
    params:
        rname = "seurat",
        sample = "{sample}",
        outdir = join(workpath, "seurat/{sample}"),
        data = join(workpath, "{sample}/outs/filtered_feature_bc_matrix/"),
        rawdata = join(workpath, "{sample}/outs/raw_feature_bc_matrix/"),
        seurat = join("workflow", "scripts", "seurat_adt.R"),
        rlibs = config["r_libs"]["ext"],
        optional = seurat_optional_params
    envmodules:
        "R/4.1"
    shell:
        """
        # Add external packages to R .libPath
        export R_LIBS_USER='{params.rlibs}'
        R --no-save --args {params.outdir} {params.data} {params.rawdata} {params.sample} {genome} {params.optional} < {params.seurat} > {log}
        """

rule seurat_rmd_report:
    input:
        join(workpath, "seurat", "{sample}", "seur_cite_cluster.rds")
    output:
        html = join(workpath, "seurat", "{sample}", "{sample}_seurat.html")
    params:
        rname = "seurat_rmd_report",
        sample = "{sample}",
        outdir = join(workpath, "seurat/{sample}"),
        seurat = join("workflow", "scripts", "seurat_adt_plot.Rmd"),
        html = join(workpath, "seurat", "{sample}", "{sample}_seurat.html")
    envmodules:
        "R/4.1"
    shell:
        """
        R -e "rmarkdown::render('{params.seurat}', params=list(workdir = '{params.outdir}', sample='{params.sample}'), output_file = '{params.html}')"
        """

rule vcf_reorder:
    input:
        vcf = vcf,
        web = expand(join(workpath, f"{sample}", "outs", "web_summary.html"), sample=lib_samples)
    output:
        vcf = join(workpath, 'demuxlet', 'vcf', 'output.vcf')
    params:
        rname = "vcf_reorder",
        reorder = join("workflow", "scripts", "reorderVCF.py"),
        bam = expand(join(workpath, "{sample}", "outs", "possorted_genome_bam.bam"), sample=lib_samples)[0]
    envmodules:
        config["tools"]["python3"]
    shell:
        """
        python3 {params.reorder} -v {input.vcf} -b {params.bam} -o {output.vcf}
        """

rule vcf_filter_blacklist:
    input:
        vcf = rules.vcf_reorder.output.vcf
    output:
        vcf = temp(join(workpath, 'demuxlet', 'vcf', 'output.filteredblacklist.recode.vcf'))
    params:
        rname = "vcf_filter_blacklist",
        blacklist = config["references"][genome]["blacklist"],
        outname = join(workpath, 'demuxlet', 'vcf', 'output.filteredblacklist')
    envmodules:
        "vcftools"
    shell:
        """
        vcftools --vcf {input.vcf} --exclude-bed {params.blacklist} --recode --out {params.outname}
        """

rule vcf_filter_quality:
    input:
        vcf = rules.vcf_filter_blacklist.output.vcf
    output:
        vcf = join(workpath, 'demuxlet', 'vcf', 'output.strict.filtered.recode.vcf')
    params:
        rname = "vcf_filter_quality",
        outname = join(workpath, 'demuxlet', 'vcf', 'output.strict.filtered')
    envmodules:
        "vcftools"
    shell:
        """
        vcftools --vcf {input.vcf} --minGQ 10 --max-missing 1.0 --recode --out {params.outname}
        """

rule demuxlet_patient_list:
    input:
        config = demuxlet_input
    output:
        patient_lists = expand(join(workpath, 'demuxlet', 'output', '{sample}', 'patient_list'), sample=lib_samples)
    params:
        rname = "demuxlet_patient_list",
        script = join("workflow", "scripts", "create_demuxlet_patient_list.py")
    envmodules:
        config["tools"]["python3"]
    shell:
        """
        python3 {params.script} {input.config}
        """

rule demuxlet_barcode:
    input:
        join(workpath, "{sample}", "outs", "web_summary.html")
    output:
        barcode = join(workpath, 'demuxlet', 'output', '{sample}', 'barcodes.tsv')
    params:
        rname = "demuxlet_barcode",
        barcodezip = join(workpath, "{sample}", "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    shell:
        """
        gunzip -c {params.barcodezip} > {output.barcode}
        """

rule run_demuxlet:
    input:
        demuxlet_input,
        vcf = rules.vcf_filter_quality.output.vcf,
        barcode = rules.demuxlet_barcode.output.barcode
    output:
        join(workpath, 'demuxlet', 'output', '{sample}', '{sample}.best'),
        join(workpath, 'demuxlet', 'output', '{sample}', '{sample}.single'),
        join(workpath, 'demuxlet', 'output', '{sample}', '{sample}.sing2')
    params:
        rname = "demuxlet",
        out = join(workpath, 'demuxlet', 'output', "{sample}", "{sample}"),
        bam = join(workpath, "{sample}", "outs", "possorted_genome_bam.bam"),
        flag = demuxlet_flag
    # envmodules: config["tools"]["python3"]
    container: config["images"]["cite_base"]
    shell:
        """
        demuxlet --group-list {input.barcode} --field GT --sam {params.bam} --vcf {input.vcf} --out {params.out} {params.flag} --alpha 0 --alpha 0.5
        """

rule seurat_aggregate:
    input:
        rds = expand(join(workpath, "seurat", "{sample}", "seur_cite_cluster.rds"), sample=lib_samples)
    output:
        rds = join(workpath, "seurat", "SeuratAggregate", "multimode.integrated.rds")
    log:
        join("seurat", "SeuratAggregate", "seurat.log")
    params:
        rname = "seurat_aggregate",
        sample = "aggregate",
        outdir = join(workpath, "seurat", "SeuratAggregate"),
        seurat = join("workflow", "scripts", "seurat_adt_aggregate.R"),
        sample_ref = lib_samples[0]
    envmodules:
        "R/4.1"
    shell:
        """
        R --no-save --args {params.outdir} {genome} {params.sample_ref} {input.rds} < {params.seurat} > {log}
        """

rule seurat_aggregate_rmd_report:
    input:
        join(workpath, "seurat", "SeuratAggregate", "multimode.integrated.rds")
    output:
        html = join(workpath, "seurat", "SeuratAggregate", "SeuratAggregate_seurat.html")
    params:
        rname = "seurat_aggregate_rmd_report",
        sample = "Aggregate",
        outdir = join(workpath, "seurat", "SeuratAggregate"),
        seurat = join("workflow", "scripts", "seurat_adt_aggregate_report.Rmd"),
        html = join(workpath, "seurat", "SeuratAggregate", "SeuratAggregate_seurat.html")
    envmodules:
        "R/4.1"
    shell:
        """
        R -e "rmarkdown::render('{params.seurat}', params=list(workdir = '{params.outdir}', sample='{params.sample}'), output_file = '{params.html}')"
        """

def azimuth_input(wildcards):
    if wildcards.sample == 'SeuratAggregate':
        return(rules.seurat_aggregate.output.rds)
    elif wildcards.sample in lib_samples:
        return(rules.seurat.output.rds)

rule azimuth:
    input:
        rds = azimuth_input
    output:
        rds = join(workpath, "azimuth", "{sample}", "azimuth_prediction.rds")
    log:
        join(workpath, "azimuth", "{sample}", "azimuth.log")
    params:
        rname = "azimuth",
        sample = "{sample}",
        outdir = join(workpath, "azimuth", "{sample}"),
        script = join("workflow", "scripts", "azimuth_celltyping.R"),
        azimuth_ref = config["references"][genome]["azimuth_ref"],
        rlibs = config["r_libs"]["ext"],
    envmodules:
        "R/4.1"
    shell:
        """
        # Add external packages to R .libPath
        export R_LIBS_USER='{params.rlibs}'
        R --no-save --args {params.outdir} {input.rds} {genome} pbmc "{params.azimuth_ref}" < {params.script} > {log}
        """
