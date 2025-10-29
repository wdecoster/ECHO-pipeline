# rules/common.smk

configfile: "config.yaml"

# Config Params
SAMPLES = config["samples"]
START_FROM = config.get("start_from", "pod5")
OUTPUT_DIR = config["output_dir"]
INPUT_DIR = config["input_dir"]
REFERENCE = config["reference"]
REFERENCE_NAME = config["reference_name"]
REFERENCE_TE = config["reference_TE"]
MIN_READ_QUAL = config["fastq_filtering"]["min_read_quality"]
MIN_READ_LENGTH = config["fastq_filtering"]["min_read_length"]
TR_CATALOG = config["tr_catalog"]
TE_CATALOG = config["te_catalog"]
FLANKING_LENGTH_BP = config["flanking_length_bp"]
CONSENSUS_EXTENSION = config["extension_repeat_consensus"]
HAPLOID_CHRS = config["haploid_chrs"]
TYPE_OF_TE = config["type_of_te"]
TYPE_OF_TR = config["type_of_tr"]


# Load packages
from pathlib import Path
import os

# Create directory to store slurm outputs
os.makedirs(f"{OUTPUT_DIR}/logs/slurm", exist_ok=True)

# Get the absolute path to the Snakefile's directory for launching custom scripts
WORKFLOW_ROOT = Path(workflow.basedir)
TR_LONGTR_METH_SCRIPT_PATH = WORKFLOW_ROOT / "scripts/TR-longTR-methylation_v4.sh"
REF_TE_METH_CPG_RES_SCRIPT_PATH = WORKFLOW_ROOT / "scripts/ref_TE_cpg_res.sh"
REF_TE_METH_AVERAGES_SCRIPT_PATH = WORKFLOW_ROOT / "scripts/ref_TE_avg_meth.sh"
TLDR_METH_SCRIPT_PATH = WORKFLOW_ROOT / "scripts/TLDR-methylation_v2.sh"


# Check if filtering is on or off (if both MIN_READ_QUAL and MIN_READ_LENGTH variables are 0, filtering is skipped)
DO_FILTER = (MIN_READ_QUAL > 0) or (MIN_READ_LENGTH > 0)

# Where to read *raw* FASTQs from (when START_FROM == "fastq").
# Falls back to the pipeline's default OUTPUT location if not provided.
FASTQ_DIR = config.get("fastq_dir", f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq")

def raw_fastq(wc):
    return f"{FASTQ_DIR}/{wc.sample}/{wc.sample}.fastq"

def fastq_for_pipeline(wc):
    if DO_FILTER:
        return f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{wc.sample}/{wc.sample}_filtered.fastq"
    else:
        return f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{wc.sample}/{wc.sample}.fastq"

#Define outputs files
all_inputs = []

if START_FROM == "pod5":
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/00_raw_data/basecalled/ubam/{{sample}}/{{sample}}_unaligned.bam", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/basecalling/{{sample}}/{{sample}}_summary.tsv", sample=SAMPLES),
    ])
elif START_FROM == "ubam":
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{{sample}}/{{sample}}.fastq", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/qc/fastq/pre_filtering/{{sample}}/{{sample}}_fastq_NanoPlot-report.html", sample=SAMPLES),
    ])
elif START_FROM == "fastq":
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/qc/fastq/pre_filtering/{{sample}}/{{sample}}_fastq_NanoPlot-report.html", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/01_alignment/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_sorted.bam", sample=SAMPLES)
    ])
elif START_FROM == "bam":
    # no basecalling or alignment outputs expected
    pass
else:
    print(f"ERROR: Invalid 'start_from' value in config.yaml: {START_FROM}")
    exit(1)


#Add filtered fastq output when filtering is applied
if START_FROM == "fastq" and DO_FILTER:
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/00_raw_data/basecalled/fastq/{{sample}}/{{sample}}_filtered.fastq", sample=SAMPLES)
    ])

#Add QC for filtered fastq if filtering was applied
if (START_FROM == "ubam" or START_FROM == "fastq") and DO_FILTER:
    all_inputs.extend([
        expand(f"{OUTPUT_DIR}/qc/fastq/post_filtering/{{sample}}/{{sample}}_filtered_fastq_NanoPlot-report.html", sample=SAMPLES)
    ])

#Default file outputs independent from which input file is used
all_inputs.extend([
    expand(f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_cramino_output.txt", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/qc/phased_bam/{{sample}}/{REFERENCE_NAME}/nanoplot/{{sample}}_{REFERENCE_NAME}_phased_bam_NanoPlot-report.html", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/qc/multiqc/{{sample}}/{REFERENCE_NAME}/multiqc_report.html", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/02_variant_calling/SNVs_Indels/{{sample}}/{REFERENCE_NAME}/phased_merge_output.vcf.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/02_variant_calling/SVs/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_SV_unphased.vcf.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/03_phasing/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_phased_alignment.bam", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_1.bed.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/phased/{{sample}}_{REFERENCE_NAME}_haplotype_2.bed.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/04_methylation_calling/{{sample}}/{REFERENCE_NAME}/unphased/{{sample}}_{REFERENCE_NAME}_unphased.bed.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/05_non_ref_TE_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}.table.txt", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/06_TR_calling/{{sample}}/{REFERENCE_NAME}/{{sample}}_{REFERENCE_NAME}_TRs_{TYPE_OF_TR}.vcf.gz", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/07_TR_methylation_calling/{{sample}}/{REFERENCE_NAME}/{TYPE_OF_TR}/{{sample}}_{REFERENCE_NAME}_TR_methylation_summary.tsv", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/non_ref_TE/{{sample}}_{REFERENCE_NAME}.table.pass.summary.meth.phased.txt", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/cpg_resolution/mod_phased/{{sample}}_{REFERENCE_NAME}_{TYPE_OF_TE}_upstream_pileup_1.bed", sample=SAMPLES),
    expand(f"{OUTPUT_DIR}/08_TE_methylation_calling/{{sample}}/{REFERENCE_NAME}/ref_TE/{TYPE_OF_TE}/{{sample}}_{REFERENCE_NAME}_methylation_summary.bed", sample=SAMPLES),
    REFERENCE
])

all_inputs.append(LOGFILE)

print(f"DEBUG: START_FROM is {START_FROM}")

