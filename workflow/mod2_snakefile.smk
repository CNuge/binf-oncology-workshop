#!/usr/bin/env python3

"""Pipeline title

AUTHORS:
    C. Nugent

Trying to turn the code from workshop into something reproducible and trackable.

### Usage

This script is run from the repository's home directory using the command :

    snakemake --snakefile workflow/mod2_snakefile.smk --cores 1

Or you can test run it with:

    snakemake --snakefile workflow/mod2_snakefile.smk -np
"""

from snakemake.utils import validate
import pandas as pd
import pathlib
import subprocess


def annotate_remote_file(fn: str):
    """
    Annotation of remote files for data transfer with provenance.
    Implicit annotation of s3 and http files based on prefix, while a
    local annotation is applied to all other (non-remote) files

    Revised version for snakemake 8+
    see: https://snakemake.github.io/snakemake-plugin-catalog/plugins/storage/s3.html
    This relies on following conda packages being in scope:
          - snakemake-storage-plugin-s3
          - snakemake-storage-plugin-http
    """
    if fn.startswith("http"):
        return storage(fn)  # noqa F821
    elif fn.startswith("s3"):
        return storage(fn)  # noqa F821
    else:
        return local(fn)  # noqa F821


# this will print the commit-ish ID to stdout for provenance
try:
    label = subprocess.check_output(["git", "describe", "--always"], encoding="UTF-8").strip()
    print(f"Precision Oncology Workshop workflow {label}")
except subprocess.CalledProcessError:
    print("Precision Oncology Workshop workflow, version not detected!")


# configure shell behavior for all rules
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

# reference the config file
configfile: "config/config.yaml"
#validate(config, schema="../schema/config.schema.yaml")

# TODO - could scale the samples with a manifest
# read in the manifest
samples = pd.read_csv(config["manifest"], sep="\t")
SAMPLES = samples["samples"].values

# --- Define path variables from config ---
BASE_DIR = config["base_dir"]
FASTQ_DIR = os.path.join(BASE_DIR, config["fastq_dir"])
PROCESSING_DIR = os.path.join(BASE_DIR, config["processing_dir"])
FINAL_BAM_DIR = os.path.join(BASE_DIR, config["final_bam_dir"])
LOG_DIR = os.path.join(BASE_DIR, config["log_dir"])
REF_DIR = os.path.join(BASE_DIR, config["ref_dir"])
REFERENCE = os.path.join(REF_DIR, config["reference"])
DBSNP = os.path.join(REF_DIR, config["dbsnp"])
GNOMAD = os.path.join(REF_DIR, config["gnomad"])

# --- Define container images from config ---
FASTQC_CONTAINER = config["containers"]["fastqc"]
BWA_CONTAINER = config["containers"]["bwa"]
SAMTOOLS_CONTAINER = config["containers"]["samtools"]
GATK_CONTAINER = config["containers"]["gatk"]



TARGETS = []

if config["s3_output_targets"]["report_bucket"]:
    TARGETS.append(
        annotate_remote_file(config["s3_output_targets"]["report_bucket"])
        + "reportable_findings_report.html"
    )
if config["s3_output_targets"]["results_bucket"]:
    TARGETS.append("results/sync/sync.results.done.flag")

rule all:
    input:
        # the first rule should define the default target files
        TARGETS,
        # Final recalibrated BAM for each sample
        expand(os.path.join(FINAL_BAM_DIR, "{sample}.recalibrated.bam"), sample=SAMPLES),
        # Initial FastQC reports
        expand(os.path.join(PROCESSING_DIR, "{sample}_test_1_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(PROCESSING_DIR, "{sample}_test_2_fastqc.html"), sample=SAMPLES)
        # Final annotated VCF
        os.path.join(VCF_DIR, "mutect2_tumour_normal.filtered.regions_of_interest.vep.txt")

# put other rules here!

# --- Helper functions for Docker paths ---
# GATK tools use '/gatk/CourseData' as the mount point
def gatk_path(path):
    return path.replace(BASE_DIR, "/gatk/CourseData")

# Other tools will use '/data' as the mount point
def data_path(path):
    return path.replace(BASE_DIR, "/data")


rule run_fastqc:
    input:
        r1=os.path.join(FASTQ_DIR, "{sample}_test_1.fastq.gz"),
        r2=os.path.join(FASTQ_DIR, "{sample}_test_2.fastq.gz")
    output:
        html1=os.path.join(PROCESSING_DIR, "{sample}_test_1_fastqc.html"),
        html2=os.path.join(PROCESSING_DIR, "{sample}_test_2_fastqc.html"),
        zip1=os.path.join(PROCESSING_DIR, "{sample}_test_1_fastqc.zip"), # Add zip to outputs
        zip2=os.path.join(PROCESSING_DIR, "{sample}_test_2_fastqc.zip")  # Add zip to outputs
    params:
        docker_mount=f"-v {BASE_DIR}:/data",
        out_dir=data_path(PROCESSING_DIR),
    shell:
        """
        docker run --rm {params.docker_mount} {FASTQC_CONTAINER} \
            fastqc -o {params.out_dir} {input.r1} {input.r2} 
        """

# -------------------------------------------------------------------
# RULE 2: BWA MEM
#
# Align paired-end reads to the reference genome.
# -------------------------------------------------------------------
rule bwa_mem:
    input:
        r1=os.path.join(FASTQ_DIR, "{sample}_test_1.fastq.gz"),
        r2=os.path.join(FASTQ_DIR, "{sample}_test_2.fastq.gz"),
        ref=REFERENCE
    output:
        sam=temp(os.path.join(PROCESSING_DIR, "{sample}.sam")) # Mark SAM as temporary
    params:
        docker_mount=f"-v {BASE_DIR}:/data",
        in_ref=data_path(REFERENCE),
    shell:
        """
        docker run --rm {params.docker_mount} {BWA_CONTAINER} \
            bwa mem {input.ref} {input.r1} {input.r2} -o {output.sam}
        """

# -------------------------------------------------------------------
# RULE 3: Samtools Sort
#
# Convert SAM to sorted, indexed BAM.
# -------------------------------------------------------------------
rule samtools_sort:
    input:
        sam=os.path.join(PROCESSING_DIR, "{sample}.sam")
    output:
        bam=os.path.join(PROCESSING_DIR, "{sample}.bam"),
        idx=os.path.join(PROCESSING_DIR, "{sample}.bam.bai")
    threads: 4
    params:
        docker_mount=f"-v {BASE_DIR}:/data",
    shell:
        """
        docker run --rm {params.docker_mount} {SAMTOOLS_CONTAINER} \
            samtools sort -@ {threads} --write-index -O bam -o {output.bam} {input.sam} 
        """

# -------------------------------------------------------------------
# RULE 4: GATK AddOrReplaceReadGroups
#
# Add read group information to the BAM file.
# -------------------------------------------------------------------
rule gatk_add_read_groups:
    input:
        bam=os.path.join(PROCESSING_DIR, "{sample}.bam"),
        idx=os.path.join(PROCESSING_DIR, "{sample}.bam.bai") # GATK needs index
    output:
        bam=temp(os.path.join(PROCESSING_DIR, "{sample}.rg.bam")) # Mark as temporary
    params:
        docker_mount=f"-v {BASE_DIR}:/gatk/CourseData",
        rgid=lambda w: f"{w.sample}_lane1",
        rglb=lambda w: f"{w.sample}_lane1",
        rgsm=lambda w: f"{w.sample}"
    shell:
        """

        docker run --rm {params.docker_mount} {GATK_CONTAINER} gatk \
            AddOrReplaceReadGroups -I {input.bam} -O {output.bam} \
            --RGID {params.rgid} --RGLB {params.rglb} --RGSM {params.rgsm} \
            --RGPL ILLUMINA --RGPU Illumina 
        """

# -------------------------------------------------------------------
# RULE 5: GATK SetNmMdAndUqTags
#
# Fix NM, MD, and UQ tags.
# -------------------------------------------------------------------
rule gatk_set_nm_md_uq_tags:
    input:
        bam=os.path.join(PROCESSING_DIR, "{sample}.rg.bam")
    output:
        bam=temp(os.path.join(PROCESSING_DIR, "{sample}.fixed.bam")) # Mark as temporary
    params:
        docker_mount=f"-v {BASE_DIR}:/gatk/CourseData",
        ref=gatk_path(REFERENCE)
    shell:
        """
        docker run --rm {params.docker_mount} {GATK_CONTAINER} gatk \
            SetNmMdAndUqTags -I {input.bam} -O {output.bam} -R {params.ref} 
        """

# -------------------------------------------------------------------
# RULE 6: GATK MarkDuplicates
#
# Mark duplicate reads.
# -------------------------------------------------------------------
rule gatk_mark_duplicates:
    input:
        bam=os.path.join(PROCESSING_DIR, "{sample}.fixed.bam")
    output:
        bam=os.path.join(PROCESSING_DIR, "{sample}.dup_marked.bam"),
        bai=os.path.join(PROCESSING_DIR, "{sample}.dup_marked.bai"), # GATK creates index
        metrics=os.path.join(PROCESSING_DIR, "{sample}.duplicate_metrics.txt")
    params:
        docker_mount=f"-v {BASE_DIR}:/gatk/CourseData",
    shell:
        """
        docker run --rm {params.docker_mount} {GATK_CONTAINER} gatk \
            MarkDuplicates -I {input.bam} -O {output.bam} \
            -M {output.metrics} --CREATE_INDEX true 
        """

# -------------------------------------------------------------------
# RULE 7: GATK BaseRecalibrator
#
# Generate the BQSR table.
# -------------------------------------------------------------------
rule gatk_base_recalibrator:
    input:
        bam=os.path.join(PROCESSING_DIR, "{sample}.dup_marked.bam"),
        bai=os.path.join(PROCESSING_DIR, "{sample}.dup_marked.bai"),
        ref=REFERENCE,
        dbsnp=DBSNP,
        gnomad=GNOMAD
    output:
        table=os.path.join(PROCESSING_DIR, "{sample}.bqsr_recal.table")
    params:
        docker_mount=f"-v {BASE_DIR}:/gatk/CourseData",
        ref=gatk_path(REFERENCE),
        dbsnp=gatk_path(DBSNP),
        gnomad=gatk_path(GNOMAD),
    shell:
        """
        docker run --rm {params.docker_mount} {GATK_CONTAINER} gatk \
            BaseRecalibrator -I {input.bam} -R {params.ref} \
            --known-sites {params.dbsnp} --known-sites {params.gnomad} \
            -O {output.table}
        """

# -------------------------------------------------------------------
# RULE 8: GATK ApplyBQSR
#
# Apply the BQSR table to create the final recalibrated BAM.
# -------------------------------------------------------------------
rule gatk_apply_bqsr:
    input:
        bam=os.path.join(PROCESSING_DIR, "{sample}.dup_marked.bam"), # Input is dup_marked BAM
        bai=os.path.join(PROCESSING_DIR, "{sample}.dup_marked.bai"),
        table=os.path.join(PROCESSING_DIR, "{sample}.bqsr_recal.table"),
        ref=REFERENCE
    output:
        bam=os.path.join(FINAL_BAM_DIR, "{sample}.recalibrated.bam"),
        bai=os.path.join(FINAL_BAM_DIR, "{sample}.recalibrated.bai") # GATK creates index
    params:
        docker_mount=f"-v {BASE_DIR}:/gatk/CourseData",
        ref=gatk_path(REFERENCE),
    shell:
        """
        docker run --rm {params.docker_mount} {GATK_CONTAINER} gatk \
            ApplyBQSR -R {params.ref} -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output.bam} --CREATE_INDEX true 
        """



# -------------------------------------------------------------------
# RULE 9: GATK Mutect2
#
# Call somatic variants using a tumour and matched normal BAM.
# Note: This rule does not use {sample} wildcard as it requires
#       both 'normal' and 'tumour' inputs explicitly.
# -------------------------------------------------------------------
rule gatk_mutect2:
    input:
        tumour_bam=os.path.join(FINAL_BAM_DIR, "tumour.recalibrated.bam"),
        tumour_bai=os.path.join(FINAL_BAM_DIR, "tumour.recalibrated.bai"),
        normal_bam=os.path.join(FINAL_BAM_DIR, "normal.recalibrated.bam"),
        normal_bai=os.path.join(FINAL_BAM_DIR, "normal.recalibrated.bai"),
        ref=REFERENCE,
        gnomad=GNOMAD
    output:
        vcf=os.path.join(VCF_DIR, "mutect2_tumour_normal.regions_of_interest.vcf.gz"),
        stats=os.path.join(VCF_DIR, "mutect2_tumour_normal.regions_of_interest.vcf.gz.stats")
    log:
        os.path.join(LOG_DIR, "gatk_mutect2", "mutect2.log")
    params:
        docker_mount=f"-v {BASE_DIR}:/gatk/CourseData",
        ref=gatk_path(REFERENCE),
        gnomad=gatk_path(GNOMAD),
    shell:
        """
        docker run --rm {params.docker_mount} {GATK_CONTAINER} gatk \
            --java-options "-Xmx12G" Mutect2 \
            -R {params.ref} \
            -I {params.input.tumour_bam)} \
            -I {params.input.normal_bam)} \
            -normal normal \
            --germline-resource {params.gnomad} \
            -O {output.vcf}
        """

# -------------------------------------------------------------------
# RULE 10: GATK FilterMutectCalls
#
# Filter the raw VCF from Mutect2.
# -------------------------------------------------------------------
rule gatk_filter_mutect_calls:
    input:
        vcf=os.path.join(VCF_DIR, "mutect2_tumour_normal.regions_of_interest.vcf.gz"),
        stats=os.path.join(VCF_DIR, "mutect2_tumour_normal.regions_of_interest.vcf.gz.stats"),
        ref=REFERENCE
    output:
        vcf=os.path.join(VCF_DIR, "mutect2_tumour_normal.filtered.regions_of_interest.vcf.gz")
    log:
        os.path.join(LOG_DIR, "gatk_filter_mutect_calls", "filter.log")
    params:
        docker_mount=f"-v {BASE_DIR}:/gatk/CourseData",
        ref=gatk_path(REFERENCE),
    shell:
        """
        docker run --rm {params.docker_mount} {GATK_CONTAINER} gatk \
            --java-options "-Xmx12G" FilterMutectCalls \
            -R {params.ref} \
            -V {input.vcf} \
            -O {output.vcf}
        """

# -------------------------------------------------------------------
# RULE 11: VEP Annotation
#
# Annotate the filtered VCF using Ensembl VEP via Apptainer/Singularity.
# -------------------------------------------------------------------
rule vep_annotate:
    input:
        vcf=os.path.join(VCF_DIR, "mutect2_tumour_normal.filtered.regions_of_interest.vcf.gz")
    output:
        txt=os.path.join(VCF_DIR, "mutect2_tumour_normal.filtered.regions_of_interest.vep.txt")
    log:
        os.path.join(LOG_DIR, "vep_annotate", "vep.log")
    params:
        sif=config["vep_sif"],
        vcf_mount=f"-B {VCF_DIR}:/data",
        cache_mount=f"-B {config['vep_cache_dir']}:/cache",
        cache_dir="/cache",
    shell:
        """
        mkdir -p $(dirname {log})
        apptainer exec {params.vcf_mount} {params.cache_mount} {params.sif} vep \
            -i {input.vcf} \
            -o {output.txt} \
            --cache \
            --dir_cache {params.cache_dir} \
            --species homo_sapiens \
            --assembly GRCh38 \
            --everything &> {log}
        """


if config["s3_output_targets"]["report_bucket"]:

    rule export_report_to_s3:
        input:
            "results/report/reportable_findings_report.html",
        output:
            annotate_remote_file(
                config["s3_output_targets"]["report_bucket"] + "reportable_findings_report.html"
            ),
        shell:
            "cp {input} {output}"


if config["s3_output_targets"]["results_bucket"]:

    rule export_pipeline_results_to_s3:
        input:
            "results/report/reportable_findings_report.html",
        output:
            "results/sync/sync.results.done.flag",
        params:
            bucket=config["s3_output_targets"]["results_bucket"].rstrip("/"),
        shell:
            "aws s3 sync results/ {params.bucket}/results/ --only-show-errors &> {output} && "
            "aws s3 sync workflow/ {params.bucket}/workflow/ --only-show-errors &>> {output} && "
            "aws s3 sync config/ {params.bucket}/config/ --only-show-errors &>> {output} && "
            "aws s3 sync resources/ {params.bucket}/resources/ --only-show-errors &>> {output} && "
            "aws s3 sync schema/ {params.bucket}/schema/ --only-show-errors &>> {output} && "
            "aws s3 cp environment.yaml {params.bucket}/environment.yaml &>> {output} && "
            "aws s3 cp README.md {params.bucket}/README.md &>> {output}"
