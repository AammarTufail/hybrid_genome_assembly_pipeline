# Antimicrobial Peptide (AMP) Detection Pipeline based on Machine Learning

## Table of Contents

- [Antimicrobial Peptide (AMP) Detection Pipeline based on Machine Learning](#antimicrobial-peptide-amp-detection-pipeline-based-on-machine-learning)
  - [Table of Contents](#table-of-contents)
  - [1. amPEPpy](#1-ampeppy)
    - [1.1 Installation](#11-installation)
    - [1.2 Training Databases](#12-training-databases)
    - [1.3 Model Training](#13-model-training)
    - [1.4 AMP Prediction](#14-amp-prediction)
    - [1.5 Results Filtering](#15-results-filtering)
    - [1.6 Output Files](#16-output-files)
  - [2. Macrel](#2-macrel)
    - [2.1 Installation](#21-installation)
    - [2.2 Step 1: Local Contig-based AMP Prediction](#22-step-1-local-contig-based-amp-prediction)
    - [2.3 Step 2: AMPSphere Exact Query](#23-step-2-ampsphere-exact-query)
    - [2.4 Step 3: AMPSphere HMMER Query](#24-step-3-ampsphere-hmmer-query)
    - [2.5 Step 4: sORF Extraction](#25-step-4-sorf-extraction)
    - [2.6 Output Files](#26-output-files)

---

## 1. amPEPpy
amPEPpy is a machine learning-based tool for predicting antimicrobial peptides (AMPs) from protein sequences. It uses a random forest classifier trained on a combined dataset of AMP and non-AMP sequences to identify potential AMPs in genomic data.
> Link for amPEPpy github repository: [https://github.com/tlawrence3/amPEPpy](https://github.com/tlawrence3/amPEPpy)

### 1.1 Installation

1. Clone the amPEPpy repository:

```bash
git clone https://github.com/tlawrence3/amPEPpy.git
cd amPEPpy
```

2. We recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html) to install amPEPpy. You can create a new conda environment and install the required dependencies with the following commands:

```bash
conda create -n amPEP python=3.8 pandas numpy biopython scikit-learn
conda activate amPEP
```

Alternatively, install using pip:

```bash
pip install .
```

3. Install amPEPpy and verify the installation:

```bash
python setup.py install
ampep -h
```

### 1.2 Training Databases

The model was trained using a combined positive dataset from three AMP databases:

| Database | Description |
|----------|-------------|
| `M_model_train_AMP_sequence.numbered.fasta` | Original amPEPpy AMP training sequences |
| `ADAPTABLE.fasta` | ADAPTABLE database (Ramos-Martín et al., 2019) — 40,000+ nonredundant sequences |
| `AMPSphere_v.2022-03.faa` | AMPSphere (Santos-Júnior et al., 2024) — 863,498 non-redundant peptide sequences |

All AMP sequences were combined into a single file: `00_training_data_all.fasta`

The negative dataset used: `M_model_train_nonAMP_sequence.numbered.fasta` contains 4,774 non-AMP sequences from the original amPEPpy training set.

### 1.3 Model Training

We trained the random forest classifier on the combined dataset:

```bash
ampep train \
    -p training_data/00_training_data_all.fasta \
    -n training_data/M_model_train_nonAMP_sequence.numbered.fasta \
    --seed 2012
```

This produces the `amPEP.model` file. 
> The pre-trained model is also available in the `pretrained_models/` directory. If you want to use the pre-trained model directly, you can skip the training step and use the provided `amPEP.model` for predictions.

**Model accuracy:**

| Metric | Value |
|--------|-------|
| Out-of-bag accuracy | 0.9958 |
| Out-of-bag balanced accuracy | 0.9894 |

### 1.4 AMP Prediction

Run predictions on Bakta-annotated protein FASTA files (`.faa`) for all genomes:

```bash
mkdir results
for i in genomes_bakta/*.faa; do
    filename=$(basename "$i" .faa)
    ampep predict -m amPEP.model -i "$i" -o "results/${filename}.tsv" --seed 2012
done
```

**Input:** Bakta-annotated `.faa` files located in `genomes_bakta/`
> You can also use Prokka or any other annotation tool to generate protein FASTA files for prediction, as long as they are in the correct format.

**Output:** Per-genome `.tsv` result files in `results/`

### 1.5 Results Filtering

Results are filtered and concatenated using the Python notebook [code.ipynb](./amPEPpy/code.ipynb). The filtering criteria:

- `probability_AMP >= 0.99`
- `probability_AMP <= 1`
- `predicted == 'AMP'`

```python
filtered_df = df[(df['probability_AMP'] >= 0.99) & 
                 (df['probability_AMP'] <= 1) & 
                 (df['predicted'] == 'AMP')]
```

Sequence information is retrieved from the original FASTA files using [extract_sequence.ipynb](./amPEPpy/extract_sequence.ipynb).

### 1.6 Output Files

```
amPEPpy/
├── amPEP.model                          # Trained random forest model
├── pretrained_models/
│   └── amPEP.model                      # Pre-trained model backup
├── genomes_bakta/
│   └── *.faa                            # Input protein FASTA files (18 genomes)
├── results/
│   ├── *.tsv                            # Per-genome prediction results (18 files)
│   └── filtered_concatenated_data.csv   # Filtered and combined AMP predictions
├── training_data/
│   └── 00_training_data_all.fasta       # Combined AMP training database
├── code.ipynb                           # Filtering and concatenation notebook
└── extract_sequence.ipynb               # Sequence retrieval notebook
```

---

## 2. Macrel

### 2.1 Installation

Install Macrel using mamba/conda:

```bash
mamba create --name env_macrel -c bioconda macrel
mamba activate env_macrel
macrel --help
```

### 2.2 Step 1: Local Contig-based AMP Prediction

Script: [`01_macrel_run.sh`](./macrel/01_macrel_run.sh)

Macrel scans assembled contig nucleotide files (`.fna`) to identify small ORFs and classify them as AMPs:

```bash
input_dir="bakta_fna_files"
master_output_dir="01_macrel_local_results"
mkdir -p "$master_output_dir"

for fasta_file in $input_dir/*.fna; do
    base_name=$(basename "$fasta_file" .fna)
    output_dir="$master_output_dir/${base_name}_results"

    macrel contigs \
        --tag "$base_name" \
        --fasta "$fasta_file" \
        --output "$output_dir" \
        --threads 32 \
        --keep-fasta-headers \
        --log-file "$output_dir/macrel_contigs.log"

    gunzip "$output_dir/${base_name}.percontigs.gz"
    gunzip "$output_dir/${base_name}.prediction.gz"
done
```

**Input:** Bakta-annotated `.fna` files in `bakta_fna_files/`

**Output:** Per-genome results in `01_macrel_local_results/` containing:
- `<genome>.prediction` — AMP classification results
- `<genome>.percontigs` — per-contig AMP count summary

### 2.3 Step 2: AMPSphere Exact Query

Script: [`02_macrel_run_amsphere.sh`](./macrel/02_macrel_run_amsphere.sh)

Queries predicted sORFs against the AMPSphere database using exact sequence matching:

```bash
input_dir="smorfs_faa_files"
master_output_dir="02_macrel_ampsphere_results"
mkdir -p "$master_output_dir"

for faa_file in $input_dir/*.smorfs.faa; do
    base_name=$(basename "$faa_file" .smorfs.faa)
    output_dir="$master_output_dir/${base_name}_ampsphere"

    macrel query-ampsphere \
        --fasta "$faa_file" \
        --tag "$base_name" \
        --output "$output_dir" \
        --log-file "$output_dir/macrel_query-ampsphere.log"

    gunzip "$output_dir/${base_name}.ampsphere_exact.tsv.gz"
done
```

**Input:** sORF `.faa` files in `smorfs_faa_files/`

**Output:** Per-genome results in `02_all_ampsphere_results/` containing:
- `<genome>.ampsphere_exact.tsv` — exact matches against AMPSphere

### 2.4 Step 3: AMPSphere HMMER Query

Script: [`03_macrel_run_amsphere_hmmer.sh`](./macrel/03_macrel_run_amsphere_hmmer.sh)

Queries predicted sORFs against the AMPSphere database using HMMER-based homology search:

```bash
input_dir="smorfs_faa_files"
master_output_dir="03__ampsphere_hmmer_results"
mkdir -p "$master_output_dir"

for faa_file in $input_dir/*.smorfs.faa; do
    base_name=$(basename "$faa_file" .smorfs.faa)
    output_dir="$master_output_dir/${base_name}_ampsphere"

    macrel query-ampsphere \
        --fasta "$faa_file" \
        --query-mode=hmmer \
        --tag "$base_name" \
        --output "$output_dir" \
        --log-file "$output_dir/macrel_query-ampsphere-hmmer.log"

    gunzip "$output_dir/${base_name}.ampsphere_hmmer.tsv.gz"
done
```

**Input:** sORF `.faa` files in `smorfs_faa_files/`

**Output:** Per-genome results in `03_all_ampsphere_hmmer_results/` containing:
- `<genome>.ampsphere_hmmer.tsv` — HMMER-based AMPSphere matches

### 2.5 Step 4: sORF Extraction

Script: [`04_macrel_run_sorfs.sh`](./macrel/04_macrel_run_sorfs.sh)

Extracts small ORFs (sORFs) from contig nucleotide files without AMP classification:

```bash
input_dir="bakta_fna_files"
master_output_dir="04_macrel_only_sORFs_results"
mkdir -p "$master_output_dir"

for fasta_file in $input_dir/*.fna; do
    base_name=$(basename "$fasta_file" .fna)
    output_dir="$master_output_dir/${base_name}_results"

    macrel get-smorfs \
        --tag "$base_name" \
        --fasta "$fasta_file" \
        --output "$output_dir" \
        --threads 32 \
        --keep-fasta-headers \
        --log-file "$output_dir/macrel_get-smorfs.log"
done
```

**Input:** Bakta-annotated `.fna` files in `bakta_fna_files/`

**Output:** Per-genome sORF FASTA files in `04_macrel_only_sORFs_results/`

### 2.6 Output Files

```
macrel/
├── 01_macrel_local_results/             # Contig-based AMP predictions
│   └── <genome>_results/
│       ├── <genome>.prediction          # AMP classification results
│       └── <genome>.percontigs          # Per-contig summary
├── 02_all_ampsphere_results/            # AMPSphere exact matches
│   └── <genome>_ampsphere/
│       └── <genome>.ampsphere_exact.tsv
├── 03_all_ampsphere_hmmer_results/      # AMPSphere HMMER matches
│   └── <genome>_ampsphere/
│       └── <genome>.ampsphere_hmmer.tsv
├── 04_macrel_only_sORFs_results/        # Extracted sORFs
│   └── <genome>_results/
│       └── <genome>.smorfs.faa
├── bakta_fna_files/                     # Input nucleotide FASTA files
├── smorfs_faa_files/                    # Input sORF protein FASTA files
├── 01_macrel_run.sh                     # Step 1 script
├── 02_macrel_run_amsphere.sh            # Step 2 script
├── 03_macrel_run_amsphere_hmmer.sh      # Step 3 script
└── 04_macrel_run_sorfs.sh              # Step 4 script
```

---

