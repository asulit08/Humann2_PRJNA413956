## =============================================================================
## WORKFLOW PROJECT: Humann2 on PRJNA413956
import pandas as pd
import glob, os, os.path, datetime, sys, csv
from os.path import join, abspath
from snakemake.utils import validate, min_version

## If SettingwithCopyWarning is ever generated, raise error instead of warning; makes code more strict
pd.set_option("mode.chained_assignment", "raise")

## =============================================================================
## set minimum snakemake version #####
## =============================================================================
min_version("5.4.0")


## =============================================================================
## SETUP
## =============================================================================

## DEFAULT PATHS
#DIR_ENVS    = "/home/arielle/conda_yaml/20200214_Humann.yaml"

## LOAD VARIABLES FROM CONFIGFILE
## submit on command-line via --configfile
if config=={}:
    sys.stderr.write("Please submit config-file with "--configfile <file>". Exit.\n")
    sys.exit(1)

## Validate configfile with yaml-schema

## define global Singularity image for reproducibility
## USE: "--use-singularity" to run all jobs in container
## without need to install any tools


## DATABASE PATHS


## Setup result dirs
DIR_BASE       = config["resultdir"]
DIR_LOGS       = join(DIR_BASE, "logs")
DIR_BENCHMARKS = join(DIR_BASE, "benchmarks")
DIR_RES        = join(DIR_BASE, "results")

#star_wrap = config["star_wrap"]

## =============================================================================
## SAMPLES
## =============================================================================

SAMPLES = pd.read_csv(abspath(config["samples"]), sep="\t", keep_default_na=False, na_values="").set_index("sample", drop=False)
##validate(SAMPLES, schema=join(DIR_SCHEMAS, "samples.schema.yaml"))

## reading samplename from samplesheet
sys.stderr.write("Reading samples from samplesheet: '{}' ...\n".format(config["samples"]))

## test if sample configuration is correct
for fname in SAMPLES["sample"]:

    fname1 = pd.notnull(SAMPLES.loc[fname, "file_mate1"])  # True or False
    fname2 = pd.notnull(SAMPLES.loc[fname, "file_mate2"]) # True or False
    fname3 = pd.notnull(SAMPLES.loc[fname, "file_singles"])  # True or False
    
    if not fname1 and not fname2:
        if not fname3:
            sys.stderr.write("No input files indicated for sample {}. Input files are necessary. Exit\n".format(fname))
            sys.exit(1)

    if fname1:
        if not os.path.isfile(SAMPLES.loc[fname, "file_mate1"]):
            sys.stderr.write("File '{}: file_mate1' from samplesheet can not be found. Make sure the file exists. Exit\n".format(fname))
            sys.exit(1)
        if not fname2:
            sys.stderr.write("For {} PE, please specify path to forward and reverse reads in columns 1 and 2, respectively\n".format(fname))
            sys.exit(1)
        else:
            if not os.path.isfile(SAMPLES.loc[fname, "file_mate2"]):
                sys.stderr.write("File '{}: file_mate2' from samplesheet can not be found. Make sure the file exists. Exit\n".format(fname))
                sys.exit(1)
    
    ## NOTE: fname1 and fname2 MUST always be present together if fname1 is true, see above and this checks if fname2 is present as well
    # if fname2 is the only 1 present, an error would be given
    if fname2:
        if not fname1:
            sys.stderr.write("For {} PE, please specify path to forward and reverse reads in columns 1 and 2, respectively".format(fname))
            sys.exit(1)
    
    if fname3:
        if not os.path.isfile(SAMPLES.loc[fname, "file_singles"]):
            sys.stderr.write("File '{}: file_singles' from samplesheet can not be found. Make sure the file exists. Exit\n".format(fname))
            sys.exit(1)

## Combinations of Input Read Files Passed Checks     
sys.stderr.write("Combinations of Input Read Files Passed Checks\n")

NUM_SAMPLES = len(SAMPLES["sample"])
sys.stderr.write("{} samples to process\n".format(NUM_SAMPLES))
sys.stderr.flush()


## =============================================================================
## WORKFLOW SETUP
## =============================================================================


## =============================================================================
## SETUP INPUT TARGETS
## =============================================================================

## Directories?


## =============================================================================
## FUNCTIONS
## =============================================================================

def get_singles_raw(wildcards):
    return SAMPLES.loc[wildcards.sample, "file_singles"]


def get_fq_pe(wildcards):
    return [SAMPLES.loc[wildcards.sample, "file_mate1"], 
            SAMPLES.loc[wildcards.sample, "file_mate2"]]


def get_data_for_mapping_pe(wildcards):
    # either raw or trimmed
    if config["trimming"]["perform"]:
        return [join(DIR_RES, "samples/trimmed/pe/{}.1.fastq.gz".format(wildcards.sample)),
                join(DIR_RES, "samples/trimmed/pe/{}.2.fastq.gz".format(wildcards.sample))]
    else:
        return get_fq_pe(wildcards) 


def get_data_for_mapping_se(wildcards):
    # either raw or trimmed
    if config["trimming"]["perform"]:
        return join(DIR_RES, "samples/trimmed/se/{}.fastq.gz".format(wildcards.sample))
    else:
        return SAMPLES.loc[wildcards.sample, "file_singles"]


def get_data(wildcards):
    fq1 = pd.notnull(SAMPLES.loc[wildcards.sample, "file_mate1"])  # True or False
    fq2 = pd.notnull(SAMPLES.loc[wildcards.sample, "file_mate2"]) # True or False
    fq3 = pd.notnull(SAMPLES.loc[wildcards.sample, "file_singles"])  # True or False
    
    if fq1 and fq2 and fq3:
        if config["mapping"]["perform"]:
            return [join(DIR_RES, "samples/map/pe/{}/Unmapped.out.mate1").format(wildcards.sample),
            join(DIR_RES, "samples/map/pe/{}/Unmapped.out.mate2").format(wildcards.sample),
            join(DIR_RES, "samples/map/se/{}/Unmapped.out.mate1").format(wildcards.sample)
            ]
        elif config["trimming"]["perform"]:
            return [join(DIR_RES, "samples/trimmed/pe/{}.1.fastq.gz".format(wildcards.sample)),
                join(DIR_RES, "samples/trimmed/pe/{}.2.fastq.gz".format(wildcards.sample)), join(DIR_RES, "samples/trimmed/se/{}.fastq.gz".format(wildcards.sample))]
        else:
            return [SAMPLES.loc[wildcards.sample, "file_mate1"], SAMPLES.loc[wildcards.sample, "file_mate2"], SAMPLES.loc[wildcards.sample, "file_singles"]]
    elif fq1 and fq2 and not fq3:
        if config["mapping"]["perform"]:
            return [join(DIR_RES, "samples/map/pe/{}/Unmapped.out.mate1").format(wildcards.sample),
            join(DIR_RES, "samples/map/pe/{}/Unmapped.out.mate2").format(wildcards.sample),
            ]
        elif config["trimming"]["perform"]:
            return [join(DIR_RES, "samples/trimmed/pe/{}.1.fastq.gz".format(wildcards.sample)),
                join(DIR_RES, "samples/trimmed/pe/{}.2.fastq.gz".format(wildcards.sample))
                ]
        else:
            return [SAMPLES.loc[wildcards.sample, "file_mate1"], SAMPLES.loc[wildcards.sample, "file_mate2"]]
    elif fq3 and not fq1 and not fq2:
        if config["mapping"]["perform"]:
            return join(DIR_RES, "samples/map/se/{}/Unmapped.out.mate1").format(wildcards.sample)
        elif config["trimming"]["perform"]:
            return join(DIR_RES, "samples/trimmed/se/{}.fastq.gz".format(wildcards.sample))
        else:
            return SAMPLES.loc[wildcards.sample, "file_singles"]
    else:
        sys.stderr.write("Wrong samplesheet.")
        sys.exit()
    

## =============================================================================
## RULES
## =============================================================================

#### TEST RUN ####
#rule all:
#    input: expand(join(DIR_RES, "{sample}.test.done"), sample=SAMPLES["sample"])
#rule test:
#    output:
#        touch("{sample}.test.done")
#### END TEST RUN ####

## Pseudo-rule to state the final targets, so that the whole runs

rule all:
    input: expand(join(DIR_RES, "Humann2_Run/{sample}.remove.done"), sample=SAMPLES["sample"])

rule trim_pe:
    input:
        fq = get_fq_pe
    output:
        trimmed1=temp(join(DIR_RES, "samples/trimmed/pe/{sample}.1.fastq.gz")), 
        trimmed2=temp(join(DIR_RES, "samples/trimmed/pe/{sample}.2.fastq.gz")),
        html=join(DIR_RES, "samples/trimmed/report/pe/{sample}.fastp.html"),
        json=join(DIR_RES, "samples/trimmed/report/pe/{sample}.fastp.json")
    log:
        join(DIR_LOGS, "samples/trimming/{sample}_pe.log")
    priority: 20
    benchmark:
        join(DIR_BENCHMARKS, "samples/trimming/{sample}_trim_pe.txt")
    params:
        extra=config["trimming"]["extra"]
    threads: 8
    conda:
        config["trimming"]["trim_env"]
    shell:
        "fastp --thread {threads} {params.extra} "
        "--in1 {input.fq[0]} --in2 {input.fq[1]} "
        "--out1 {output.trimmed1} --out2 {output.trimmed2} "
        "--html {output.html} --json {output.json} > {log} 2>&1"

rule trim_se:
    input:
        sample=get_singles_raw
    output:
        trimmed=temp(join(DIR_RES, "samples/trimmed/se/{sample}.fastq.gz")),
        html=join(DIR_RES, "samples/trimmed/report/se/{sample}.fastp.html"),
        json=join(DIR_RES, "samples/trimmed/report/se/{sample}.fastp.json")
    log:
        join(DIR_LOGS, "samples/trimming/{sample}_se.log")
    priority: 20
    benchmark:
        join(DIR_BENCHMARKS, "samples/trimming/{sample}_trim_se.txt")
    params:
        extra=config["trimming"]["extra"]
    threads: 8
    conda:
        config["trimming"]["trim_env"]
    shell:
        "fastp --thread {threads} {params.extra} --in1 {input.sample} "
        "--out1 {output.trimmed} "
        "--html {output.html} --json {output.json} > {log} 2>&1"

rule star_mapping_pe:
    input:
        sample=get_data_for_mapping_pe
    output:
        unmapped1=temp(join(DIR_RES, "samples/map/pe/{sample}/Unmapped.out.mate1")),
        unmapped2=temp(join(DIR_RES, "samples/map/pe/{sample}/Unmapped.out.mate2")),
        #sam=temp(join(DIR_RES, "samples/map/pe/{sample}/Aligned.out.sam"))
    priority: 10
    log:
        join(DIR_LOGS, "samples/star_mapping/{sample}_pe.log")
    benchmark:
        join(DIR_BENCHMARKS, "samples/star_mapping/{sample}_pe.txt")
    conda:
        config["mapping"]["map_env"]
    params:
        # path to STAR reference genome index
        index=config["mapping"]["index"],
        #annotation=config["mapping"]["annotation"],
        extra=config["mapping"]["star_extra"]
    threads: 8
    wrapper:
        config["mapping"]["star_wrap"]


rule star_mapping_se:
    input:
        sample=get_data_for_mapping_se
    output:
        unmapped=temp(join(DIR_RES, "samples/map/se/{sample}/Unmapped.out.mate1")),
        #sam=temp(join(DIR_RES, "samples/map/se/{sample}/Aligned.out.sam"))
    priority: 10
    log:
        join(DIR_LOGS, "samples/star_mapping/{sample}_se.log")
    benchmark:
        join(DIR_BENCHMARKS, "samples/star_mapping/{sample}_se.txt")
    conda:
        config["mapping"]["map_env"]
    params:
        # path to STAR reference genome index
        index=config["mapping"]["index"],
        #annotation=config["mapping"]["annotation"],
        extra=config["mapping"]["star_extra"]
    threads: 8
    wrapper:
        config["mapping"]["star_wrap"]

rule combine:
    input:
        get_data
    output:
        temp(join(DIR_RES, "samples/cat_input/{sample}.fq.gz"))
    benchmark:
        join(DIR_BENCHMARKS, "samples/{sample}.combineinput.txt")
    shell:
        "cat {input} | gzip > {output}"

rule humann:
    input: 
        join(DIR_RES, "samples/cat_input/{sample}.fq.gz")
    output:
        temp(touch(join(DIR_RES, "Humann2_Run/{sample}.humann2.done")))
    log:
        humout = join(DIR_LOGS, "Humann2_Run/{sample}.humann2.out"),
        humerr = join(DIR_LOGS, "Humann2_Run/{sample}.humann2.err")
    params:
        prescreen = config["prescreen"],
        outp = directory(join(DIR_RES, "Humann2_Run/{sample}_humann2_out/"))
    benchmark:
        join(DIR_BENCHMARKS, "humann2/{sample}.humann2Run.txt")
    conda:
        config["humann_env"]
    threads: 8
    #resources:
    #    const=1
    shell:
        "humann2 --input {input} --threads {threads} --prescreen-threshold {params.prescreen} --output {params.outp} 1> {log.humout} 2> {log.humerr}"

rule humann_regroup:
    input: 
        rules.humann.output
    output:
        temp(join(DIR_RES, "Humann2_Run/{sample}_humann2_out/regroup/{sample}_uniref90go.tsv"))
    params:
        inp = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/{sample}_genefamilies.tsv"),
        groups = "uniref90_go"
    log:
        regout = join(DIR_LOGS, "Humann2_Run/{sample}.humann2_regroup.out"),
        regerr = join(DIR_LOGS, "Humann2_Run/{sample}.humann2_regroup.err")
    benchmark:
        join(DIR_BENCHMARKS, "humann2/{sample}.humann2_regroup.txt")
    conda:
        config["humann_env"]
    shell:
        "humann2_regroup_table --input {params.inp} --groups {params.groups} --output {output} 1> {log.regout} 2> {log.regerr}"
    
rule humann_rename:
    input:
        rules.humann_regroup.output
    output:
        join(DIR_RES, "Humann2_Run/{sample}_humann2_out/regroup/{sample}_uniref90go_names.tsv")
    log:
        nameout = join(DIR_LOGS, "Humann2_Run/{sample}.humann2_rename.out"),
        nameerr = join(DIR_LOGS, "Humann2_Run/{sample}.humann2_rename.err")
    benchmark:
        join(DIR_BENCHMARKS, "humann2/{sample}.humann2_rename.txt")
    conda:
        config["humann_env"]
    shell:
        "humann2_rename_table -i {input} -n go -o {output} 1> {log.nameout} 2> {log.nameerr}"
    
rule humann_mv:
    input: 
        rules.humann_rename.output
    params:
        log = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/{sample}_humann2_temp/{sample}.log"),
        bugs = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/{sample}_humann2_temp/{sample}_metaphlan_bugs_list.tsv"),
        outmv = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/")
    output:
        temp(touch(join(DIR_RES, "Humann2_Run/{sample}.mv.done")))
    benchmark:
        join(DIR_BENCHMARKS, "humann2/{sample}.humann2mv.txt")
    shell:
        "mv {params.log} {params.bugs} {params.outmv}"

rule humann_compress:
    input:
        rules.humann_mv.output
    params:
        outcom = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/{sample}_humann2_temp.tar.gz"),
        outdir = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/"),
        in2comp = "{sample}_humann2_temp"
    output:
        temp(touch(join(DIR_RES, "Humann2_Run/{sample}.compress.done")))
    benchmark:
        join(DIR_BENCHMARKS, "humann2/{sample}.humann2compress.txt")
    shell:
        "tar -C {params.outdir} -cvzf {params.outcom} {params.in2comp}"

rule humann_rm:
    input:
        rules.humann_compress.output
    output:
        temp(touch(join(DIR_RES, "Humann2_Run/{sample}.remove.done")))
    params:
        torem = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/{sample}_humann2_temp")
    benchmark:
        join(DIR_BENCHMARKS, "humann2/{sample}.humann2rm.txt")
    shell:
        "rm -r {params.torem}"
