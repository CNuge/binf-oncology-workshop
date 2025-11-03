#!/bin/bash

set -euo pipefail

DATE=$(date +"%Y%m%d%H%M")

target=$(grep "results_bucket" config/config.yaml | cut -d":" -f2- | sed 's/"//g' | sed 's/\/$//' | sed 's/.*~.*//')


## NOTE: change the filenames below to reflect the pipeline being developed

if [[ -z "${target}" ]]; then
    snakemake --resources load=150 --profile ../sb-slurm/ &> seqbio_gwas_qc_${DATE}.log
else
    snakemake --resources load=150 --profile ../sb-slurm/ &>> seqbio_gwas_qc_${DATE}.log &&
        aws s3 cp seqbio_gwas_qc_${DATE}.log $target/seqbio_gwas_qc_${DATE}.log
fi
