import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

sample = snakemake.input.get("sample")
# fq1 = snakemake.input.get("fq1")
assert sample is not None, "input-> fq1 is a required input parameter"
fq1 = (
    [snakemake.input[0]] if isinstance(snakemake.input[0], str) else snakemake.input[0]
)

if len(snakemake.input) > 1:
    fq2 = snakemake.input[1]
    fq2 = (
        [snakemake.input[1]]
        if isinstance(snakemake.input[1], str)
        else snakemake.input[1]
    )
    assert len(fq1) == len(
        fq2
    ), "input-> equal number of files required for fq1 and fq2"
else:
    fq2 = None

input_str_fq1 = ",".join(fq1)
input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
input_str = " ".join([input_str_fq1, input_str_fq2])

if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

outprefix = os.path.dirname(snakemake.output[0]) + "/"

shell(
    "STAR "
    "{extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.params.index} "
    "--readFilesIn {input_str} "
    "{readcmd} "
    "--outFileNamePrefix {outprefix} "
    "--outReadsUnmapped Fastx "
    "--outSAMmode None "
    "--outSAMtype None "
    "--outStd Log "
    "{log}"
)

# "--outSAMmode None " # if no sam is desired
#     "--quantMode GeneCounts "
