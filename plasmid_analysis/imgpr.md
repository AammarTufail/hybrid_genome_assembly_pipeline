<!-- Created by Hendrik Bethge, 06.05.2026 -->

# IMG/PR Plasmid Identification using BLASTn

## Table of Contents

- [IMG/PR Plasmid Identification using BLASTn](#imgpr-plasmid-identification-using-blastn)
  - [Table of Contents](#table-of-contents)
  - [1. Overview](#1-overview)
  - [2. BLAST Installation](#2-blast-installation)
  - [3. IMG/PR Database Preparation](#3-imgpr-database-preparation)
  - [4. Plasmid Identification using BLASTn](#4-plasmid-identification-using-blastn)
  - [5. Metadata Integration](#5-metadata-integration)
  - [6. Output Files](#6-output-files)

---

## 1. Overview

This workflow identifies plasmids by comparing predicted plasmid sequences against the IMG/PR plasmid database using **BLASTn**.

For every query plasmid sequence, the pipeline searches for highly similar plasmids in IMG/PR and retrieves associated metadata from the IMG/PR database.

The workflow consists of:

- Creating a local BLAST database from IMG/PR sequences
- Running BLASTn against the IMG/PR database
- Filtering hits using sequence identity and query coverage thresholds
- Integrating IMG/PR metadata with BLAST results

The resulting table can be used to identify plasmids in IMG/PR that are similar to the query plasmid sequences and retrieve their associated metadata.

> Link for IMG/PR database: https://genome.jgi.doe.gov/portal/IMG_PR

> Link for Blastn: https://www.ncbi.nlm.nih.gov/books/NBK279690/

---

## 2. BLAST Installation

We recommend using Micromamba.

Create a dedicated environment and install BLAST:

```bash
micromamba create -n blastn -c conda-forge -c bioconda blast -y

micromamba activate blastn

blastn -version
```

---

## 3. IMG/PR Database Preparation

Download the IMG/PR plasmid database from https://genome.jgi.doe.gov/portal/IMG_PR/IMG_PR.home.html and create a local BLAST database.

```bash
mkdir -p imgpr_db
cd imgpr_db

unzip IMG_PR_download.zip

gunzip IMG_VR_2023-08-08_1/IMGPR_nucl.fna.gz

makeblastdb \
    -in IMG_VR_2023-08-08_1/IMGPR_nucl.fna \
    -dbtype nucl \
    -out IMGPR_blastdb
```
---

## 4. Plasmid Identification using BLASTn

The workflow compares plasmid sequences against the IMG/PR plasmid database.

### Input

```text
# Blast database
IMGPR_blastdb.* 

# Query: Plasmid sequences
plasmids.fasta

# IMG/PR metadata
IMG_VR_2023-08-08_1/IMGPR_plasmid_data.tsv
```

### Run Analysis

```bash
blastn \
    -query plasmids.fasta \
    -db IMGPR_blastdb \
    -out imgpr_results.tsv \
    -outfmt 6 \
    -num_threads 8 \
    -perc_identity 85 \
    -evalue 1e-10 \
    -qcov_hsp_perc 50
```

### BLAST Parameters

| Parameter | Value | Description |
|------------|---------|-------------|
| `perc_identity` | 85 | Minimum nucleotide identity (%) |
| `qcov_hsp_perc` | 50 | Minimum query coverage (%) |
| `evalue` | 1e-10 | Maximum accepted E-value |
| `outfmt` | 6 | Tabular output format |

The analysis was performed using a minimum nucleotide identity of 85% and a minimum query coverage of 50% to allow detection of more distant plasmid homologs.

### Output

```text
imgpr_results.tsv
```
Each row corresponds to a BLAST hit against an IMG/PR plasmid sequence.

---

## 5. Metadata Integration

IMG/PR metadata are added to the BLAST results using the IMG/PR plasmid identifiers.

```bash
# extract IMG/PR plasmid ID from BLAST subject column
awk -F"\t" '{split($2,a,"|"); print $0 "\t" a[1]}' \
    imgpr_results.tsv > blast_with_id.tsv

# sort BLAST results
sort -k13,13 blast_with_id.tsv > blast_sorted.tsv

# prepare metadata table
cut -f1,2- IMG_VR_2023-08-08_1/IMGPR_plasmid_data.tsv \
    > IMG_VR_2023-08-08_1/metadata_prepared.tsv

# sort metadata
sort -k1,1 IMG_VR_2023-08-08_1/metadata_prepared.tsv \
    > IMG_VR_2023-08-08_1/metadata_sorted.tsv

# join metadata and BLAST results
join -t $'\t' \
    -1 13 \
    -2 1 \
    blast_sorted.tsv \
    IMG_VR_2023-08-08_1/metadata_sorted.tsv \
    > blast_with_metadata.tsv
```
"blast_with_metadata.tsv" represents the final result of the workflow and contains BLAST alignments together with the corresponding IMG/PR metadata.

### Output

```text
blast_with_metadata.tsv
```

This file contains BLAST hits together with IMG/PR metadata.

---

## 6. Output Files

```text
IMGPR/
├──imgpr_db/
│   ├── IMG_VR_2023-08-08_1/
│   │   ├── IMGPR_nucl.fna
│   │   ├── README.md
│   │   ├── IMGPR_plasmid_data.tsv
│   │   ├── metadata_prepared.tsv
│   │   └── metadata_sorted.tsv
│   │
│   └── IMGPR_blastdb.*                 # BLAST database files
│
├── plasmids.fasta                   # Plasmid fasta file
│
├── imgpr_results.tsv                # Raw BLAST results
├── blast_with_id.tsv                # BLAST results with IMG/PR IDs
├── blast_sorted.tsv                 # Sorted BLAST results
├── blast_with_metadata.tsv          # Final annotated BLAST results
│
└── lower_parameters_ident85_qcov50/
    ├── imgpr_results.tsv
    ├── blast_with_id.tsv
    ├── blast_sorted.tsv
    └── blast_with_metadata.tsv
```

---
