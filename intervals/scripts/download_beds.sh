#!/bin/bash

gsutil cp \
    "gs://broad-dsde-methods-palantir-data/AoU_gene_reports/hg38_maybe/acmg59_allofus_19dec2019.GRC38.wGenes.NEW.bed" .

gsutil cp \
    "gs://broad-dsde-methods-palantir-data/AoU_gene_reports/hg38_maybe/GRCh38_lowmappabilityall.bed.gz" .

gsutil cp \
    "gs://broad-dsde-methods-palantir-data/AoU_gene_reports/hg38_maybe/GRCh38_gc85_slop50.bed.gz" .

gsutil cp \
    "gs://broad-dsde-methods-palantir-data/AoU_gene_reports/hg38_maybe/GRCh38_gclt25_merged.bed" .

gsutil -u "$(gcloud config get-value project)" cp \
    "gs://broad-dsde-methods-hydro-gen-truth-data-public/NIST/GIAB/v4/giab_highconf_wgs_calling_regions_hg38_intersection.bed" .
