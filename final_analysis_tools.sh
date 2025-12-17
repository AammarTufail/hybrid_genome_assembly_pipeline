#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=1-23:55:00
#SBATCH --job-name=analysis_tools
#SBATCH --output=analysis_tools.out
#SBATCH --error=analysis_tools.err
#SBATCH --partition=base

# -----------------------------------------------------------------------------
# ------------------------------ Info -----------------------------------------
# -----------------------------------------------------------------------------
# Pipeline to install and run the genome analysis tools
    # Please change: 
        # home_dir (line 42)
        # db_home_dir (line 45)
        # Micromamba activation (line 28-30)
    # Place the annotated assemblies and reads in the respective directories before running the pipeline (line 50-51)
    
#--------------Tools:
# TXSScan (https://github.com/macsy-models/TXSScan)
    # Identification and annotation of bacterial secretion systems
# dbCAN3 (https://github.com/linnabrown/run_dbcan)
    # identification and annotation of CAZymes and PULs 
# IMG/PR (https://img-dev.jgi.doe.gov/cgi-bin/plasmid/main.cgi)
    # Plamid database used to identify potential plasmids assembled by Plassembler
    # Blast performed with Blastn (https://github.com/JacobLondon/Blastn)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Load modules
module load micromamba/1.4.2
export MAMBA_ROOT_PREFIX=$HOME/.micromamba
eval "$(micromamba shell hook --shell=bash)"

# define the maximum amount of threads available
threads=$(nproc)

# ##-------------------------create the folder structure--------------------
# Home directory
home_dir="/path/to/your/assembly_pipeline_directory"
      #change to your desired home directory for the assembly pipeline
db_home_dir="/path/to/your/databases_directory"
      #if changed the gtdb database download directory needs to be updated in line: 142
mkdir -p $db_home_dir
cd $home_dir || exit 1

# --- Input_directory ---
short_input_dir="$home_dir/00_raw_reads/raw_short_reads"
long_input_dir="$home_dir/00_raw_reads/raw_long_reads"
genomes_fasta="$home_dir/all_assemblies/fasta"
genomes_gff="$home_dir/all_bakta"
    #place all raw reads in the correct folder

# --- Output Directories ---
# TXSScan output
txsscan_out="$home_dir/txsscan"
combined_file="$txsscan_out/combined_results.tsv" 
combined_file_gff="$txsscan_out/combined_results_annotated.tsv"
mkdir -p $txsscan_out
# dbCAN3 Output
dbcan3_out="$home_dir/dbcan3"
mkdir -p $dbcan3_out
# IMG/PR blast output
imgpr_out="$home_dir/imgpr"
mkdir -p $imgpr_out
#---------------------------------------------------------------------------------------
# ##-----------------------------Install & Run the analysis Tools-----------------------
#---------------------------------------------------------------------------------------

######################################################################################
## ------------------------------------ TXSScan --------------------------------------
# ---------------install 
micromamba create -n txsscan python=3.13.3 -y
micromamba activate txsscan 
micromamba install bioconda::macsyfinder=2.1.4 -y

# Download TXSScan database
mkdir -p $db_home_dir//txsscan
macsydata install TXSscan -t $db_home_dir/txsscan # version 1.1.3
micromamba deactivate

# ---------------- Run TXSScan
micromamba activate txsscan0
cd $txsscan_out || exit 1

# Run txsscan
for i in $genomes_fasta/*.fasta; do
  name=$(basename "$i" .faa)

  out_dir="$txsscan_out/${name}"
  mkdir -p "$out_dir"
  
  macsyfinder --models TXSscan all \
  --models-dir $db_home_dir/txsscan/ \
  --sequence-db $i \
  --db-type unordered \
  -o "$out_dir" 
done

#---- combine into one results file
# Initialize output file with header (assuming all_systems.tsv has a header)
  # Add "sample" column to the header
first_file=true # to grab the header only once

for file in $txsscan_out/*/all_systems.tsv; do
   sample=$(basename $(dirname "$file"))
   
       # Skip empty or non-existent files
   [ ! -s "$file" ] && continue
   
   # Get header line (first line with at least 3 columns, excluding comments)
   header=$(grep -v '^#' "$file" | awk -F '\t' 'NF >= 3 {print; exit}')

   if [ "$first_file" = true ]; then
       # Write header + extra column name to combined file
       echo -e "${header}\tsample" > "$combined_file"
       first_file=false
   fi

   # Append data: skip comment lines and any rows matching the header
   grep -v '^#' "$file" | awk -v sample="$sample" -v header="$header" -F '\t' '
       NF >= 3 && $0 != header { print $0 "\t" sample }
   ' >> "$combined_file"
done

#------ Add annotation information + contig and position from gff3
# Print header (TXSScan header + GFF3 header)
{
    head -n 1 "$combined_file" | tr '\n' '\t'
    echo -e "gff_contig\tgff_source\tgff_type\tgff_start\tgff_end\tgff_score\tgff_strand\tgff_phase\tgff_attributes"
} > $combined_file_gff

# Process each TXSScan line (skipping header)
tail -n +2 "$combined_file" | while read -r line; do
    sample=$(echo "$line" | awk -F'\t' '{print $NF}')
    hit_id=$(echo "$line" | awk -F'\t' '{print $2}')
    gff_file="$genomes_gff/$sample/${sample}.gff3"

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
        echo -e "${line}\t${match_line}" >> "$combined_file_gff"
    else
        echo "No GFF3 match for $hit_id in $gff_file" >&2
    fi
done

micromamba deactivate

######################################################################################
# ##---------------------------------- Run dbcan3 ------------------------------------
# ----- Install dbcan3 -----
micromamba create -n dbcan3 python=3.8 -y
micromamba activate dbcan3
micromamba install bioconda::dbcan=4.1.4 -y
micromamba install matplotlib=3.4.2 -y
micromamba install -c conda-forge libnsl=2.0.0 -y 
# Blastp runs into issues with libnsl.so.1 missing, conda forge does not provide libnsl.so.1, only newer version
  # create a symlink from libnsl.so.3 to libnsl.so.1
    ln -s $HOME/.micromamba/envs/dbcan3/lib/libnsl.so.3 $HOME/.micromamba/envs/dbcan3/lib/libnsl.so.1
    export LD_LIBRARY_PATH=$HOME/.micromamba/envs/dbcan3/lib:$LD_LIBRARY_PATH
# Check installations
blastp -h
run_dbcan -h

# Download dbcan3 database
dbcan_build --cpus 8 --db-dir $db_home_dir/dbcan_db --clean
    # doesnt work properly, so need to create the blastdb manually afterwards
makeblastdb -in $db_home_dir/dbcan_db/PUL.faa -dbtype prot

micromamba deactivate dbcan3

# ----- Run dbcan3 -----
micromamba activate dbcan3
cd $dbcan3_out || exit 1

# run dbcan 
for i in $genomes_fasta/*.fasta; do
  
    # Get the base name of the file (without the path and extension)
    base=$(basename "$i" _assembly.fasta)

    outdir="$dbcan3_out/$base"
    mkdir -p $outdir
    
    run_dbcan $i prok --out_dir $outdir --db_dir $dbcan_db --gram n --cluster cluster --cgc_substrate --hmm_cpu 16 --dbcan_thread 16
done

#-----------------------------------------Create summary files-----------------------------------
cd $dbcan3_out|| exit 1

# Create summary file for cgc_standard.out
output_cgc_standard="summary_cgc_standard.tsv"
  # Remove the output file if it exists, to avoid appending to an old file
rm -f $output_cgc_standard

# Loop through all directories and find cgc_standard.out
for file in $(find . -type f -name "cgc_standard.out"); do
    # Get the sample name which is the name of the directory containing the file
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


# Create summary file for cgc_standard.out
output_substrate="summary_substrate.tsv"
  # Remove the output file if it exists, to avoid appending to an old file
rm -f $output_substrate

# Loop through all directories and find substrate.out
for file in $(find . -type f -name "substrate.out"); do
    # Get the sample name which is the name of the directory containing the file
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

# Create summary file for cgc_standard.out
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

micromamba deactivate

#######################################################################################
# ## --------------------------------- IMG/VR plasmid database Blastn ------------------------
micromamba create -n blastn blast=2.15.0 -y
micromamba activate blastn

# ----------- create blast db
mkdir $db_home_dir/imgpr_blastdb
cd $db_home_dir/imgpr_blastdb || exit 1

# Download IMG/PR plasmid database manually from: https://genome.jgi.doe.gov/portal/pages/accessDenied.jsf?state=%27anonDownload%27
    # register with ORCID ID
# Place the downloaded zip file in the imgpr_blastdb folder
unzip IMG_PR_download.zip -d IMG_PR
gunzip IMG_VR_2023-08-08_1/*.fna.gz IMG_VR_2023-08-08_1/
makeblastdb -in IMG_VR_2023-08-08_1/*.fna -dbtype nucl -out IMGPR_plasmids

# -------------------- blast
cd $imgpr_out || exit 1

# Run blastn against plasmid database
  # Place all circular plasmids sequences in a fasta file named: plassembler_circular_plasmids.fasta
  # change perc_identity, evalue and qcov_hsp_perc as needed
blastn -query plassembler_circular_plasmids.fasta -db $db_home_dir/IMGPR_plasmids -out imgpr_results.tsv -outfmt 6 -num_threads 8 \
 -perc_identity 95 -evalue 1e-10 -qcov_hsp_perc 80

# ----------- add metadata to blast results
awk -F"\t" '{split($2,a,"|"); print $0 "\t" a[1]}' imgpr_results.tsv > blast_with_id.tsv

# Keep only relevant columns and ensure tab-delimited
cut -f1,2- IMG_VR_2023-08-08_1/IMGPR_plasmid_data.tsv > IMG_VR_2023-08-08_1/metadata_prepared.tsv
 
# Sort BLAST by last column (plasmid ID)
sort -k13,13 blast_with_id.tsv > blast_sorted.tsv

# Sort metadata by plasmid ID (first column)
sort -k1,1 IMG_VR_2023-08-08_1/metadata_prepared.tsv > IMG_VR_2023-08-08_1/metadata_sorted.tsv

# join metadata with blast
join -t $'\t' -1 13 -2 1 blast_sorted.tsv IMG_VR_2023-08-08_1/metadata_sorted.tsv > blast_with_metadata.tsv


micromamba deactivate
# ##-----------------------------jobinfo + deactivate mamba and unload all packages-------------------
micromamba deactivate
module purge
jobinfo 