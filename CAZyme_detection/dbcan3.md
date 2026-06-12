<!-- Created by Hendrik Bethge, 06.05.2026 -->

# # CAZyme, CGC and PUL-based Substrate Prediction using dbCAN3

## Table of Contents

- - [CAZyme, CGC and PUL-based Substrate Prediction using dbCAN3](#cazyme-cgc-and-pul-based-substrate-prediction-using-dbcan3)
  - [Table of Contents](#table-of-contents)
  - [1. Overview](#1-overview)
  - [2. dbCAN3 Installation](#2-dbcan3-installation)
  - [3. Database Preparation](#3-database-preparation)
  - [4. dbCAN3 Analysis](#4-dbcan3-analysis)
  - [5. Create Summary Files](#5-create-summary-files)
  - [6. Output Files](#6-output-files)

---

## 1. Overview

This workflow identifies **Carbohydrate-Active Enzymes (CAZymes)** and **Carbohydrate Gene Clusters (CGCs)** using **dbCAN3**.

For each genome, dbCAN3 predicts CAZyme families and Carbohydrate Gene Clusters (CGCs). Identified CGCs are subsequently compared against the dbCAN PUL database to infer potential carbohydrate substrates.

The workflow consists of:

- Running dbCAN3 on genome FASTA files
- Identification of CAZyme families
- Detection of Carbohydrate Gene Clusters (CGCs)
- Comparison of CGCs against known Polysaccharide Utilization Loci (PULs)
- Prediction of CGC substrates
- Combination of results from multiple genomes into summary tables

The resulting tables can be used to compare CAZyme repertoires and carbohydrate utilization potential across genomes.

> Link for dbCAN3: https://github.com/linnabrown/run_dbcan

> Link for dbCAN database: https://bcb.unl.edu/dbCAN2/

> Link for CAZy databases: https://www.cazy.org

---

## 2. dbCAN3 Installation

We recommend using Micromamba.

Create a dedicated environment and install dbCAN3:

```bash
micromamba create -n dbcan3 python=3.8 -y

micromamba activate dbcan3

micromamba install bioconda::dbcan=4.1.4 -y

micromamba install matplotlib=3.4.2 -y

micromamba install -c conda-forge libnsl=2.0.0 -y
```

### Additional BLAST Compatibility Fix
Installations may fail because `blastp` expects `libnsl.so.1`, while Conda only provides newer library versions.

To fix this create a symlink from libnsl.so.3 to libnsl.so.1
```bash
ln -s \
$HOME/.micromamba/envs/dbcan3/lib/libnsl.so.3 \
$HOME/.micromamba/envs/dbcan3/lib/libnsl.so.1

export LD_LIBRARY_PATH=$HOME/.micromamba/envs/dbcan3/lib:$LD_LIBRARY_PATH
```

Verify the installation:

```bash
blastp -h

run_dbcan -h
```

---

## 3. Database Preparation

Download the dbCAN database files.

```bash
dbcan_build \
    --cpus 8 \
    --db-dir dbcan_db \
    --clean
```

The dbCAN download does not always generate the required BLAST database automatically. Therefore, a manual database creation step may be required:

```bash
makeblastdb \
    -in dbcan_db/PUL.faa \
    -dbtype prot
```

### Output

```text
dbcan_db/
├── PUL.faa
├── PUL.faa.phr
├── PUL.faa.pin
├── PUL.faa.psq
└── ...
```

---

## 4. dbCAN3 Analysis

The workflow runs dbCAN3 on each genome independently.

### Input

```text
# dbCAN database
dbcan_db/

# Genome files
genomes_fasta/
├── SampleA.fasta
├── SampleB.fasta
└── ...
```

### Run Analysis

```bash
micromamba activate dbcan3

cd dbcan3_out || exit 1

for i in genomes_fasta/*.fasta; do

    sample=$(basename "$i" _assembly.fasta)

    output_dir="dbcan3_out/${sample}"
    mkdir -p "$output_dir"

    run_dbcan \
        "$i" \
        prok \
        --out_dir "$output_dir" \
        --db_dir dbcan_db \
        --gram n \
        --cluster cluster \
        --cgc_substrate \
        --hmm_cpu 16 \
        --dbcan_thread 16

done

micromamba deactivate
```

### Output

Each genome receives its own result directory:

```text
dbcan3_out/
├── SampleA/
│   ├── overview.txt
│   ├── cgc_standard.out
│   ├── substrate.out
│   └── ...
├── SampleB/
│   ├── overview.txt
│   ├── cgc_standard.out
│   ├── substrate.out
│   └── ...
└── ...
```

---

## 5. Create Summary Files

Results from all genomes are combined into three summary tables.

### Summary of CGC Predictions

Input files:

```text
*/cgc_standard.out
```

Create summary file for cgc_standard.out
```bash
output_cgc_standard="summary_cgc_standard.tsv"
rm -f $output_cgc_standard
    # Remove the output file if it exists, to avoid appending to an old file

# Loop through all directories and find cgc_standard.out
for file in $(find . -type f -name "cgc_standard.out"); do
    # Define the sample name (name of the directory)
    sample_name=$(basename $(dirname "$file"))
    
     # If the output file doesn't have the header yet, add the header from the first file
    if [ ! -f $output_cgc_standard ]; then
        # Copy the header (first line) from the first file to the output
        head -n 1 "$file" | sed 's/^/genome_name\t/' >> "$output_cgc_standard"
    fi
    
    # Append the rest of the lines (excluding the header) to the output file
    tail -n +2 "$file" | while read -r line
    do
        # Append the sample name and the line to the summary file
        echo -e "$sample_name\t$line" >> "$output_cgc_standard"
    done
done

```
Output:

```text
summary_cgc_standard.tsv
```

This file contains all predicted carbohydrate gene clusters (CGCs) across all genomes.

---

### Summary of Substrate Predictions

Input files:

```text
*/substrate.out
```
During processing, the original CGC identifier is split into:

- `contig`
- `cgc`

to facilitate downstream analysis.

Create summary file for cgc_standard.out
```bash
output_substrate="summary_substrate.tsv"
rm -f $output_substrate
    # Remove the output file if it exists, to avoid appending to an old file

# Loop through all directories and find substrate.out
for file in $(find . -type f -name "substrate.out"); do
    # Define the sample name (name of the directory)
    sample_name=$(basename $(dirname "$file"))
    
    # If the output file doesn't have the header yet, add the modified header from the first file
    if [ ! -f $output_substrate ]; then
        # Copy the header (first line) from the first file and modify it to include the new columns "contig" and "cgc"
        head -n 1 "$file" | sed 's/^#cgcid/contig\tcgc/' | sed 's/^/genome_name\t/' >> "$output_substrate"
    fi
    
    # Append the rest of the lines (excluding the header) to the output file
    tail -n +2 "$file" | while read -r line; do
        # Split the first column (#cgcid) into two columns: "contig" and "cgc"
        contig=$(echo "$line" | cut -d'|' -f1)
        cgc=$(echo "$line" | cut -d'|' -f2)
        
        # Remove the original #cgcid column and append "contig", "cgc", and the rest of the line to the output
        rest_of_line=$(echo "$line" | cut -d'|' -f2- --complement)
        
        # Append the sample name, contig, cgc, and the rest of the line to the summary file
        echo -e "$sample_name\t$contig\t$cgc\t$rest_of_line" >> "$output_substrate"
    done
done
``` 

Output:

```text
summary_substrate.tsv
```
This file contains substrate predictions for all identified CGCs.

---

### Summary of Genome-wide CAZyme Profiles

Input files:

```text
*/overview.txt
```

Create summary file for cgc_standard.out
```bash
output_overview="summary_overview.tsv"
  # Remove the output file if it exists, to avoid appending to an old file
rm -f $output_overview

# Loop through all directories and find cgc_standard.out
for file in $(find . -type f -name "overview.txt"); do
    # Get the sample name which is the name of the directory containing the file
    sample_name=$(basename $(dirname "$file"))
    
     # If the output file doesn't have the header yet, add the header from the first file
    if [ ! -f $output_overview ]; then
        # Copy the header (first line) from the first file to the output
        head -n 1 "$file" | sed 's/^/genome_name\t/' >> "$output_overview"
    fi
    
    # Append the rest of the lines (excluding the header) to the output file
    tail -n +2 "$file" | while read -r line
    do
        # Append the sample name and the line to the summary file
        echo -e "$sample_name\t$line" >> "$output_overview"
    done
done
```

Output:

```text
summary_overview.tsv
```
This file contains genome-wide CAZyme annotation summaries for all genomes.



---

## 6. Output Files

```text
dbcan3_out/
├── summary_cgc_standard.tsv      # Combined CGC predictions
├── summary_substrate.tsv         # Combined substrate predictions
├── summary_overview.tsv          # Combined CAZyme overview
│
├── SampleA/
│   ├── overview.txt
│   ├── cgc_standard.out
│   ├── substrate.out
│   └── ...
│
├── SampleB/
│   ├── overview.txt
│   ├── cgc_standard.out
│   ├── substrate.out
│   └── ...
│
└── ...
```

---
