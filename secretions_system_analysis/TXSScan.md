<!-- Created by Hendrik Bethge, 06.05.2026 -->

# Bacterial Secretion System Detection using TXSScan

## Table of Contents

- [Bacterial Secretion System Detection using TXSScan](#bacterial-secretion-system-detection-using-txsscan)
  - [Table of Contents](#table-of-contents)
  - [1. Overview](#1-overview)
  - [2. TXSScan Installation](#2-txsscan-installation)
  - [3. TXSScan Database Preparation](#3-txsscan-database-preparation)
  - [4. Secretion System Detection](#4-secretion-system-detection)
  - [5. Combine Results](#5-combine-results)
  - [6. Add Genome Annotation Information](#6-add-genome-annotation-information)
  - [7. Output Files](#7-output-files)

---

## 1. Overview

This workflow identifies bacterial secretion systems using **TXSScan**, a component of the **MacSyFinder** framework designed for the detection of protein secretion systems in bacterial genomes.

For each genome, TXSScan searches for conserved secretion system components and classifies them into known secretion system families.

The workflow consists of:

- Running TXSScan on genome FASTA files
- Identification of secretion system components and complete systems
- Combination of results from multiple genomes
- Integration of gene annotation information from GFF3 files

The resulting tables can be used to compare secretion system distributions across genomes and link detected secretion system components to annotated genomic features.

> Link for TXSScan: https://github.com/macsy-models/TXSScan

> Link for MacSyFinder: https://github.com/gem-pasteur/macsyfinder

---

## 2. TXSScan Installation

We recommend using Micromamba.

Create a dedicated environment and install MacSyFinder:

```bash
micromamba create -n txsscan python=3.13.3 -y

micromamba activate txsscan

micromamba install bioconda::macsyfinder=2.1.4 -y

macsyfinder --help
```

---

## 3. TXSScan Database Preparation

Install the TXSScan models used by MacSyFinder.

```bash
mkdir -p txsscan_db

macsydata install TXSscan -t txsscan_db
```

---

## 4. Secretion System Detection

The workflow runs TXSScan on each genome FASTA file independently.

### Input

```text
# database
txsscan_db/
└── TXSscan/*

# Genome files
genomes_fasta/
├── SampleA.fasta
├── SampleB.fasta
└── ...
```

### Run Analysis

```bash
micromamba activate txsscan

cd txsscan_out || exit 1

for genome in genomes_fasta/*.fasta; do

    sample=$(basename "$genome" .fasta)

    output_dir="txsscan_out/${sample}"
    mkdir -p "$output_dir"

    macsyfinder \
        --models TXSscan all \
        --models-dir txsscan_db/ \
        --sequence-db "$genome" \
        --db-type unordered \
        -o "$output_dir"

done

micromamba deactivate
```

### Output

Each genome receives its own result directory:

```text
txsscan_out/
├── SampleA/
│   └── all_systems.tsv
├── SampleB/
│   └── all_systems.tsv
└── ...
```

The file `all_systems.tsv` contains all secretion systems detected within a genome.

---

## 5. Combine Results

Results from all genomes are combined into a single table and an additional sample column is added.

```bash
first_file=true

for file in txsscan_out/*/all_systems.tsv; do

    sample=$(basename "$(dirname "$file")")

    [ ! -s "$file" ] && continue

    header=$(grep -v '^#' "$file" | awk -F '\t' 'NF >= 3 {print; exit}')

    if [ "$first_file" = true ]; then
        echo -e "${header}\tsample" > combined_file.tsv
        first_file=false
    fi

    grep -v '^#' "$file" | \
    awk -v sample="$sample" -v header="$header" -F '\t' '
        NF >= 3 && $0 != header {
            print $0 "\t" sample
        }
    ' >> combined_file.tsv

done
```

### Output

```text
combined_file
``` 
This file contains all detected secretion systems from all genomes.

---

## 6. Add Genome Annotation Information

### Input
Gene annotation information is retrieved from genome GFF3 files and added to the combined TXSScan results.

```text
# GFF3 file for each genome
genomes_gff/
├── SampleA/
│   └── SampleA.gff3
├── SampleB/
│   └── SampleB.gff3
└── ...

# Combined TXSScan results
txsscan_out/combined_file.tsv                   
```


### Run Annotation Integration

```bash
{
    head -n 1 combined_file.tsv | tr '\n' '\t'
    echo -e "gff_contig\tgff_source\tgff_type\tgff_start\tgff_end\tgff_score\tgff_strand\tgff_phase\tgff_attributes"
} > combined_file_gff.tsv

tail -n +2 combined_file.tsv | while read -r line; do

    sample=$(echo "$line" | awk -F'\t' '{print $NF}')
    hit_id=$(echo "$line" | awk -F'\t' '{print $2}')

    gff_file="genomes_gff/$sample/${sample}.gff3"

    if [[ ! -f "$gff_file" ]]; then
        echo "Missing GFF3: $gff_file" >&2
        continue
    fi

    match_line=$(awk -v id="$hit_id" -F'\t' '
    $0 ~ /^[^#]/ && NF == 9 && $9 ~ ("ID=" id "(;|$)") {
        print $1, $2, $3, $4, $5, $6, $7, $8, $9;
        exit;
    }
    ' OFS='\t' "$gff_file")

    if [[ -n "$match_line" ]]; then
        echo -e "${line}\t${match_line}" >> combined_file_gff.tsv
    fi

done
```

### Output

```text
combined_file_gff.tsv
```

This file contains TXSScan results together with genomic annotation information from the corresponding GFF3 files.

---

## 7. Output Files

```text
txsscan_out/
├── combined_file.tsv                    # Combined TXSScan results
├── combined_file_gff.tsv                # Combined results with GFF3 annotations
│
├── SampleA/
│   ├── all_systems.tsv
│   └── ...
│
├── SampleB/
│   ├── all_systems.tsv
│   └── ...
│
└── ...
```

---
