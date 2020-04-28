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
DIR_ENVS    = "/home/arielle/conda_yaml/20200214_Humann.yaml"

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

#def get_fastq1(wildcards):
#    return SAMPLES.loc[wildcards.sample, "file_mate1"]


#def get_fastq2(wildcards):
#    return SAMPLES.loc[wildcards.sample, "file_mate2"]


#def get_singles(wildcards):
#    return SAMPLES.loc[wildcards.sample, "file_singles"]


def get_data(wildcards):
    fq1 = pd.notnull(SAMPLES.loc[wildcards.sample, "file_mate1"])  # True or False
    fq2 = pd.notnull(SAMPLES.loc[wildcards.sample, "file_mate2"]) # True or False
    fq3 = pd.notnull(SAMPLES.loc[wildcards.sample, "file_singles"])  # True or False
    
    if fq1 and fq2 and fq3:
        data = [SAMPLES.loc[wildcards.sample, "file_mate1"], SAMPLES.loc[wildcards.sample, "file_mate2"], SAMPLES.loc[wildcards.sample, "file_singles"]]
    elif fq1 and fq2 and not fq3:
        data = [SAMPLES.loc[wildcards.sample, "file_mate1"], SAMPLES.loc[wildcards.sample, "file_mate2"]]
    elif fq3 and not fq1 and not fq2:
        data = [SAMPLES.loc[wildcards.sample, "file_singles"]]
    else:
        sys.stderr.write("Wrong samplesheet.")
        sys.exit()

    return data 
    

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
    input: expand("{sample}.remove.done", sample=SAMPLES["sample"])

rule combine:
    input:
        get_data
    output:
        join(DIR_RES, "cat_input/{sample}.fq.gz")
    benchmark:
        join(DIR_BENCHMARKS, "{sample}.combineinput.txt")
    shell:
        "cat {input} | gzip > {output}"

rule humann:
    input: 
        join(DIR_RES, "cat_input/{sample}.fq.gz")
    output:
        directory(join(DIR_RES, "Humann2_Run/{sample}_humann2_out/"))
    log:
        humout = join(DIR_LOGS, "Humann2_Run/{sample}.humann2.out"),
        humerr = join(DIR_LOGS, "Humann2_Run/{sample}.humann2.err")
    benchmark:
        join(DIR_BENCHMARKS, "{sample}.humann2Run.txt")
    conda:
        DIR_ENVS
    shell:
        "humann2 --input {input} --output {output} 1> {log.humout} 2> {log.humerr}"
    
rule humann_mv:
    input: 
        rules.humann.output
    params:
        log = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/{sample}_humann2_temp/{sample}.log"),
        bugs = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/{sample}_humann2_temp/{sample}_metaphlan_bugs_list.tsv"),
        outmv = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/")
    output:
        temp(touch("{sample}.mv.done"))
    benchmark:
        join(DIR_BENCHMARKS, "{sample}.humann2mv.txt")
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
        temp(touch("{sample}.compress.done"))
    benchmark:
        join(DIR_BENCHMARKS, "{sample}.humann2compress.txt")
    shell:
        "tar -C {params.outdir} -cvzf {params.outcom} {params.in2comp}"

rule humann_rm:
    input:
        rules.humann_compress.output
    output:
        temp(touch("{sample}.remove.done"))
    params:
        torem = join(DIR_RES, "Humann2_Run/{sample}_humann2_out/{sample}_humann2_temp")
    benchmark:
        join(DIR_BENCHMARKS, "{sample}.humann2rm.txt")
    shell:
        "rm -r {params.torem}"
