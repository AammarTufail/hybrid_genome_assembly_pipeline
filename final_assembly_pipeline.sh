#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=1-23:55:00
#SBATCH --job-name=assembly_pipeline
#SBATCH --output=assembly_pipeline.out
#SBATCH --error=assembly_pipeline.err
#SBATCH --partition=base

# -----------------------------------------------------------------------------
# ------------------------------ Info -----------------------------------------
# -----------------------------------------------------------------------------
# Installation, database download and folder structure for the assembly pipeline
    # Please change: 
        # home_dir (line 65)
        # db_home_dir (line 67)
        # Micromamba activation (line 58-61)
    # Place the raw reads in the short_input_dir (lie 84) and long_input_dir (line 85), or change the path
    
# GTDB-tk could run into pplacer problems with insufficient memory.
    # in this case allocate more (at least 64 GB), or run gtdb-tk multiple times with less samples each

#--------------Tools:
# FastQC (https://github.com/s-andrews/FastQC)
    # Quality control Illumina short reads
# Fastp (https://github.com/s-andrews/FastQC)
    # Processing of Illumina short reads
# NanoPlot (https://github.com/wdecoster/NanoPlot)
    # Quality control of ONT long reads
# Filtlong (https://github.com/rrwick/Filtlong)
    # Processing of ONT long reads
# Canu (https://github.com/marbl/canu)
    # Error correction of ONT long reads
# Unicycler (https://github.com/rrwick/Unicycler
    # Hybrid genome assembly from short and long reads
# Quast (https://github.com/ablab/quast)
    # Quality assesment of the assembled genomes
# CheckM (https://github.com/Ecogenomics/CheckM)
    # Quality assesment of the assembled genomes 
        # Completeness and Contamination estimation based on marker genes
# CheckM2 (https://github.com/chklovski/CheckM2)
    # Quality assesment of the assembled genomes
        # Completeness and Contamination estimation based on machine learning
# Prokka (https://github.com/tseemann/prokka)
    # Annotation of the assembled genomes
# Bakta (https://github.com/oschwengers/bakta)
    # Annotation of the assembled genomes (more comprehensive)
# GTDB-tk (https://github.com/Ecogenomics/GTDBTk)
    # Taxonomic classification of the assembled genomes
# Plassembler (https://github.com/gbouras13/plassembler)
    # Plasmid assembly from short and long reads

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Load modules
module load micromamba/1.4.2
export MAMBA_ROOT_PREFIX=$HOME/.micromamba
eval "$(micromamba shell hook --shell=bash)"

# ##-------------------------create the folder structure--------------------
# Home directory
home_dir="/path/to/your/assembly_pipeline_directory"
      #change to your desired home directory for the assembly pipeline
db_home_dir="/path/to/your/databases_directory"
      #if changed the gtdb database download directory needs to be updated in line: 142
mkdir -p $db_home_dir
cd $home_dir || exit 1

# Databases
checkm_db_dir="$db_home_dir/checkm_db"
checkm2_db_dir="$db_home_dir/checkm2_db"
# quast_db_dir="$db_home_dir/quast_db"
    # not downloadable anymore until Quast recieves a new update for the Busco database
bakta_db_dir="$db_home_dir/bakta_db"
gtdb_db_dir="$db_home_dir/gtdb_db" 
      #if anything is changed, you need to update the path in line 225
plassembler_db_dir="$db_home_dir/plassembler_db"

# Assembly pipeline
# Input_directory
short_input_dir="$home_dir/00_raw_reads/raw_short_reads"
long_input_dir="$home_dir/00_raw_reads/raw_long_reads"
mkdir -p $short_input_dir $long_input_dir 
    #place all raw reads in the correct folder
# 01-Short read trimming + qc
fastqc_out_raw="$home_dir/01_cleaning_short_reads/01_short_reads_qc_raw"
fastp_out="$home_dir/01_cleaning_short_reads/02_cleaned_short_reads"
fastqc_out_cleaned="$home_dir/01_cleaning_short_reads/03_short_reads_qc_cleaned"
mkdir -p $fastqc_out_raw $fastp_out $fastqc_out_cleaned
# 02-Long read trimming + qc
nanoplots_out_raw="$home_dir/02_cleaning_long_reads/01_long_reads_qc_raw"
filtlong_out="$home_dir/02_cleaning_long_reads/02_cleaned_long_reads/cleaned_fastq"
filtlong_out_fasta="$home_dir/02_cleaning_long_reads/02_cleaned_long_reads/cleaned_fasta"
nanoplots_out_cleaned="$home_dir/02_cleaning_long_reads/03_long_reads_qc_cleaned"
canu_out="$home_dir/02_cleaning_long_reads/04_corrected_long_reads"
canu_corrected_reads="$canu_out/all_corrected_reads"
nanoplots_out_canu="$home_dir/02_cleaning_long_reads/05_long_reads_qc_corrected"
final_long_reads="$home_dir/02_cleaning_long_reads/06_final_long_reads"
nanoplots_out_final="$home_dir/02_cleaning_long_reads/07_long_reads_qc_final"
mkdir -p $nanoplots_out_raw $filtlong_out $filtlong_out_fasta $nanoplots_out_cleaned $canu_out $canu_corrected_reads $nanoplots_out_canu $final_long_reads $nanoplots_out_final
# 03-Assemblies
unicycler_out="$home_dir/03_assembly/01_unicycler_assemblies"
all_assemblies_fasta="$unicycler_out/all_assemblies_fasta"
all_assemblies_gfa="$unicycler_out/all_assemblies_gfa"
mkdir -p $unicycler_out $all_assemblies_fasta $all_assemblies_gfa
# 04-Assembly qc  
quast_out="$home_dir/04_assembly_qc/quast"
checkm_out="$home_dir/04_assembly_qc/checkm"
checkm2_out="$home_dir/04_assembly_qc/checkm2"
mkdir -p $quast_out $checkm_out $checkm2_out 
# 05-Annotation
prokka_out="$home_dir/05_annotated_genomes/prokka"
all_annotated_gff="$prokka_out/all_gff"
all_annotated_fsa="$prokka_out/all_fsa"
mkdir -p $prokka_out $all_annotated_gff $all_annotated_fsa
bakta_out="$home_dir/05_annotated_genomes/bakta"
bakta_gff="$bakta_out/all_gff"
bakta_fsa="$bakta_out/all_fsa"
mkdir -p $bakta_out $bakta_gff $bakta_fsa
# 06-Taxonomic classification
gtdb_out="$home_dir/06_taxonomic_classification"
mkdir -p $gtdb_out 
# 07-Plassembler
plassembler_out="$home_dir/07_plasmid_assembly"
mkdir -p $plassembler_out 

# status report file name
report=status_report_assembly.txt
  # Create the file or continue the last created one
  if [ ! -e "$report" ]; then
      echo -e "Starting assembly pipeline!" > "$report"
        # If it doesn't exist, create the file and add the starting line
  else
      printf "\n" >> $report
      echo -e "\nStarting assembly pipeline!" >> "$report"
         # If it exists, append the starting line
  fi
printf "\n" >> $report

#---------------------------------------------------------------------------------------
# ##-----Set the amount of Ram and CPUs available per task running parallel-------------
#---------------------------------------------------------------------------------------

# define the maximum amount of threads available
threads=$(nproc)

# Total number of input files
num_files_short=$(ls $short_input_dir/* | wc -l)
num_files_long=$(ls $long_input_dir/* | wc -l)

# Number of threads per job on short reads, rounded down
threads_per_job_short=$(echo "scale=0; $threads / $num_files_short" | bc)
  if [ "$threads_per_job_short" -lt 1 ]; then
    threads_per_job_short=1
  fi
# Number of threads per job on long reads/assembly, rounded down
threads_per_job_long=$(echo "scale=0; $threads / $num_files_long" | bc)
  if [ "$threads_per_job_long" -lt 1 ]; then
    threads_per_job_long=1
  fi

# create the status update file "$report"
cd $home_dir
echo "Number of threads available= $threads" >> $report
echo "Number of short reads= $num_files_short" >> $report
echo "Number of long reads= $num_files_long" >> $report
echo "Threads available for parallel jobs on short reads= $threads_per_job_short" >> $report
echo "Threads available for parallel jobs on long reads/assembly= $threads_per_job_long" >> $report
printf "\n" >> $report

#---------------------------------------------------------------------------------------
# ##---------------------------------Run the assembly pipeline--------------------------
#---------------------------------------------------------------------------------------

# #-------------------------------------01 short read trimming & qc-----------------
micromamba activate 01_short_reads_qc
#-----------------01-01 fastqc_before fastp

# Initialize a counter
counter=0

# Iterate over each .fastq.gz file in the raw reads directory
for file in $short_input_dir/*.fastq.gz
do
  # Run FastQC in the background and increment the counter
  fastqc -o ${fastqc_out_raw} $file -t $threads_per_job_short &
  counter=$((counter+1))

  # If counter has reached max_processes, wait until all background processes have finished before continuing
  if [[ $((counter%num_files_short)) -eq 0 ]]; then
    wait
    counter=0
  fi
done
# Wait for all background processes to finish before continuing
wait

# Run multiqc to aggregate the results
cd $fastqc_out_raw || exit 1
multiqc .

# print the status update and output into $report
cd $home_dir
echo -e "01-01 Finished quality control of the raw short reads: Fastqc. The output files are:" >> $report
for i in $fastqc_out_raw/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

#-------------------01-02 fastp for trimming and cleaning

cd $fastp_out

# Iterate over each pair of .fastq.gz files in the directory
for file1 in $short_input_dir/*_1.fastq.gz
do
  # Extract base of filename for output
  base=$(basename "$file1" _1.fastq.gz)

  # Define paired file
  file2=${file1/_1.fastq.gz/_2.fastq.gz}
  
  # Define the output file names
  output1=$fastp_out/$(basename ${file1})
  output2=$fastp_out/$(basename ${file2})
  
  # Run fastp in the background and increment the counter
  fastp -i $file1 -I $file2 -o $output1 -O $output2 -j $fastp_out/${base}.json -h $fastp_out/${base}.html -q 30 --dedup --thread $threads_per_job_short &
  counter=$((counter+1))

  # If counter has reached max_instances, wait until all background processes have finished before continuing
  if [[ $((counter%num_files_short)) -eq 0 ]]; then
    wait
    counter=0
  fi
done
# Wait for all background processes to finish before continuing
wait

# print the status update and output into $report
cd $home_dir
echo -e "01-02 Finished short read cleaning: Fastp. The output files are:" >> $report
for i in $fastp_out/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

#---------------------01-03 fastqc_after cleaning

# Initialize a counter
counter=0

# Iterate over each .fastq.gz file in the raw reads directory
for file in $fastp_out/*.fastq.gz
do
  # Run FastQC in the background and increment the counter
  fastqc -o $fastqc_out_cleaned $file -t $threads_per_job_short &
  counter=$((counter+1))

  # If counter has reached max_processes, wait until all background processes have finished before continuing
  if [[ $((counter%num_files_short)) -eq 0 ]]; then
    wait
    counter=0  
  fi
done
# Wait for all background processes to finish before continuing
wait

# Run multiqc to aggregate the results
cd $fastqc_out_cleaned || exit 1
multiqc .

# print the status update and output into $report
cd $home_dir
echo -e "01-03 Finished quality control of the trimmed short reads: Fastqc. The output files are:" >> $report
for i in $fastqc_out_cleaned/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report


micromamba deactivate

# ##-------------------------------02 long read trimming + qc--------------------------------------------------
micromamba activate 02_long_reads_qc
##---------------------02-01 long reads quality control before filtlong

# ## Nanoplot
# Initialize job counter
counter=0

# Loop over all the .fastq.gz files in the input directory
for file in $long_input_dir/*.fastq.gz
do
    # Get the base name of the file (without the path and extension)
    base=$(basename "$file" .fastq.gz)

    # Run NanoPlot on the file in the background
    NanoPlot --fastq $file -o $nanoplots_out_raw/$base -t $threads_per_job_long --plots kde --format png --N50 --dpi 300 --store --raw --tsv_stats --info_in_report &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
      wait
      counter=0
    fi
done
# Wait for all background jobs to finish
wait

# ## Nanocomp
# Create one array off all fastq.gz files and one for the sample names.
fastq_list=(); names_list=()
for file in $long_input_dir/*.fastq.gz; do
    sample_name=$(basename "$file" .fastq.gz)
    fastq_list+=("$file")
    names_list+=("$sample_name")
done
# Join the array elements with spaces to form the list of fastq.gz files
fastq_list="${fastq_list[*]}"
names_list="${names_list[*]}"

# run nanocomp on the list
NanoComp -t $threads --fastq $fastq_list --names $names_list --outdir $nanoplots_out_raw --raw

# print the status update and output into $report
cd $home_dir
echo -e "02-01 Finished quality control of the raw long reads: Nanoplot & Nanocomp. The output files are:" >> $report
for i in $nanoplots_out_raw/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

##---------------------02-02 long read filtering and trimming

# Initialize job counter
counter=0

# Loop over all the .fastq.gz files in the input directory
for file in $long_input_dir/*.fastq.gz
do
    # Get the base name of the file (without the path and extension)
    base=$(basename "$file" .fastq.gz)

    # Run Filtlong on the file
    filtlong --min_length 1000 --keep_percent 90 $file | gzip > $filtlong_out/${base}_filtlong.fastq.gz &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait


for i in $filtlong_out/*_filtlong.fastq.gz
do
    # Get the base name of the file (without the path and extension)
    base=$(basename "$i" _filtlong.fastq.gz)

    seqtk seq -A $i | gzip > $filtlong_out_fasta/${base}_filtlong.fasta.gz
done


# print the status update and output into $report
cd $home_dir
echo -e "02-02 Finished long read trimming: Filtlong. The output files are:" >> $report
for i in $filtlong_out/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

##---------------------02-03 long reads quality control after cleaning

## Nanoplot
# Initialize job counter
counter=0

# Loop over all the .fastq.gz files in the input directory
for file in $filtlong_out/*_filtlong.fastq.gz
do
    # Get the base name of the file (without the path and extension)
    base=$(basename "$file" _filtlong.fastq.gz)

    # Run NanoPlot on the file in the background
    NanoPlot --fastq $file -o $nanoplots_out_cleaned/$base -t $threads_per_job_long --plots kde --format png --N50 --dpi 300 --store --raw --tsv_stats --info_in_report &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# ## Nanocomp
# Create one array off all fastq.gz files and one for the sample names.
fastq_list=(); names_list=()
for file in $filtlong_out/*.fastq.gz; do
    sample_name=$(basename "$file" .fastq.gz)
    fastq_list+=("$file")
    names_list+=("$sample_name")
done
# Join the array elements with spaces to form the list of fastq.gz files
fastq_list="${fastq_list[*]}"
names_list="${names_list[*]}"

# run nanocomp on the list
NanoComp -t $threads --fastq $fastq_list --names $names_list --outdir $nanoplots_out_cleaned --raw

# print the status update and output into $report
cd $home_dir
echo -e "02-03 Finished quality control of the trimmed long reads: Nanoplot & Nanocomp. The output files are:" >> $report
for i in $nanoplots_out_cleaned/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

##---------------------02-04 long reads error correction

for i in $filtlong_out/*_filtlong.fastq.gz; do
    # Get the base name of the file (without the path and extension)
    base=$(basename "$i" _filtlong.fastq.gz)

    # Create an output directory for each sample
    canu_out_sample="$canu_out/$base"
    mkdir -p $canu_out_sample

    # Run Canu. Can't be run as a parallel loop, because is is submitting scripts that could interefere with each other
    canu -p $base -d $canu_out_sample genomeSize=4.8m maxInputCoverage=100 -nanopore $i stopAfter=trimming useGrid=false maxThreads=$threads  
done

# Copy all corrected fasta files files into one dir
find $canu_out/ -type f -name "*.correctedReads.fasta.gz" -exec cp {} $canu_corrected_reads \;

# print the status update and output into $report
cd $home_dir
echo -e "02-04 Finished error correction of the long reads: Canu. The output files are:" >> $report
for i in $canu_corrected_reads/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

##---------------------02-05 long reads quality control after read correction

# Initialize job counter
counter=0

# Loop over all the .fastq.gz files in the input directory
for i in $canu_corrected_reads/*.correctedReads.fasta.gz
do
    # Get the base name of the file (without the path and extension)
    base=$(basename "$i" .correctedReads.fasta.gz)

    # Run NanoPlot on the file in the background
    NanoPlot --fasta $i -o $nanoplots_out_canu/$base -t $threads_per_job_long --plots kde --format png --N50 --dpi 300 --store --raw --tsv_stats --info_in_report &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# ## Nanocomp
# Create one array off all fastq.gz files and one for the sample names.
fastq_list=(); names_list=()
for file in $canu_corrected_reads/*.correctedReads.fasta.gz; do
    sample_name=$(basename "$file" *.correctedReads.fasta.gz)
    fastq_list+=("$file")
    names_list+=("$sample_name")
done
# Join the array elements with spaces to form the list of fastq.gz files
fastq_list="${fastq_list[*]}"
names_list="${names_list[*]}"

# run nanocomp on the list
NanoComp -t $threads --fasta $fastq_list --names $names_list --outdir $nanoplots_out_canu --raw

# print the status update and output into $report
cd $home_dir
echo -e "02-06 Finished quality control of the corrected long reads: Nanoplot & Nanocomp. The output files are:" >> $report
for i in $nanoplots_out_canu/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

##---------------------02-06 copying the final long reads into one folder (use canu corrected if possible)

# create a list for the names of all long reads without canu correction
not_corrected_reads=()

for i in $filtlong_out_fasta/*_filtlong.fasta.gz; do
    
    # Get the base name of the file (without the path and extension)
    base=$(basename "$i" _filtlong.fasta.gz)

    # Define corresponding long read file (prefer canu_corrected over just filtlong)
    if [ -f "$canu_corrected_reads/${base}.correctedReads.fasta.gz" ]; then
        cp $canu_corrected_reads/${base}.correctedReads.fasta.gz $final_long_reads/${base}_final_longread.fasta.gz
    elif [ -f "$filtlong_out_fasta/${base}_filtlong.fasta.gz" ]; then
        cp $filtlong_out_fasta/${base}_filtlong.fasta.gz $final_long_reads/${base}_final_longread.fasta.gz
      # Store the names of samples without canu_corrected_long_reads in a list to add to the report.txt file
      not_corrected_reads+=("$base")
    fi
done

# print the status update and output into $report
cd $home_dir
echo -e "02-07 The final long reads are now stored in $final_long_reads:" >> $report
for i in $final_long_reads/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done

echo -e "For the following samples just the fastp cleaned long reads are used, since canu correction did not work:" >> $report
for i in "${not_corrected_reads[@]}"; do
    echo "$item" >> $report
done
printf "\n" >> $report
##---------------------02-07 long reads quality control of the final long reads

# Initialize job counter
counter=0

# Loop over all the .fastq.gz files in the input directory
for i in $final_long_reads/* ;do
    # Get the base name of the file (without the path and extension)
      base=$(basename "$i" _final_longread.fasta.gz)

    # Run NanoPlot on the file in the background
    NanoPlot --fasta $i -o $nanoplots_out_final/$base -t $threads_per_job_long --plots kde --format png --N50 --dpi 300 --store --raw --tsv_stats --info_in_report &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# ## Nanocomp
# Create one array off all fastq.gz files and one for the sample names.
fastq_list=(); names_list=()

# Loop over all the final long reads
for i in $final_long_reads/* ;do
    # Get the base name of the file 
      sample_name=$(basename "$i" _final_longread.fasta.gz)
    fastq_list+=("$i")
    names_list+=("$sample_name")
done
# Join the array elements with spaces to form the list of fastq.gz files
fastq_list="${fastq_list[*]}"
names_list="${names_list[*]}"

# run nanocomp on the list
NanoComp -t $threads --fasta $fastq_list --names $names_list --outdir $nanoplots_out_final --raw

# print the status update and output into $report
cd $home_dir
echo -e "02-06 Finished quality control of the final long reads: Nanoplot & Nanocomp. The output files are:" >> $report
for i in $nanoplots_out_final/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

micromamba deactivate
# ##-------------------------------03 hybrid assembly---------------------------------------------------
micromamba activate 03_assembly
#---------------------------------03_01 Unicycler

# Initialize job counter
counter=0

# Create an array with the base names of the short read files
base_names=($(ls $fastp_out/ | grep '_1.fastq.gz' | sed 's/_1.fastq.gz//g'))

# Loop over all the base names
for base in "${base_names[@]}"
do
    # Define the paired-end files
    file1=$fastp_out/${base}_1.fastq.gz
    file2=$fastp_out/${base}_2.fastq.gz
    long_read=$final_long_reads/${base}_final_longread.fasta.gz


    # Run Unicycler on the paired-end files and the long reads
    unicycler -1 $file1 -2 $file2 -l $long_read -o $unicycler_out/$base --min_fasta_length 1000 -t $threads_per_job_long &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

#rename the assembly.fasta file according to the sample name and copy into fasta collection dir
cd $unicycler_out || exit 1
  # Iterate over directories
for dir in */; do
  cd "$dir" || continue
  # Check if the assembly.fasta file exists
  if [ -f "assembly.fasta" ]; then
    # Get the name of the directory
    dir_name=$(basename "$dir")
    # copy the assembly.fasta file with sample name to fasta folder
    cp assembly.fasta $all_assemblies_fasta/"${dir_name}_assembly.fasta"
  fi
  # Navigate back to the parent directory
  cd ..
done

#rename the assembly.gfa file according to the sample name and copy into fasta collection dir
#used for bandage visualisation
cd $unicycler_out || exit 1
  # Iterate over directories
for dir in */; do
  cd "$dir" || continue
  # Check if the assembly.fasta file exists
  if [ -f "assembly.gfa" ]; then
    # Get the name of the directory
    dir_name=$(basename "$dir")
    # copy the assembly.fasta file with sample name (dir name) to gfa folder
    cp assembly.gfa $all_assemblies_gfa/"${dir_name}_assembly.gfa"
  fi
  # Navigate back to the parent directory
  cd ..
done

# print the status update and output into $report
cd $home_dir
echo -e "03_01 Finished hybrid assembly: Unicycler. The fasta output files are:" >> $report
for i in $all_assemblies_fasta/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
echo -e "Assemblies without canu corrected long reads:= $no_corrected_reads" >> $report
printf "\n" >> $report

# ##----------------------------------------------------------04 Assembly quality control--------------------------------------------------
# ##-----------------------------------------------------04-01 Quast
micromamba activate 04_quast  

cd $all_assemblies_fasta

# Run QUAST on the assembly files
quast.py -o $quast_out *.fasta --circos --glimmer --rna-finding --conserved-genes-finding --report-all-metrics --use-all-alignments -t $threads

# print the status update and output into $report
cd $home_dir
echo -e "04-01 Finished assembly quality control: Quast. The output files are:" >> $report
for i in $quast_out/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report
   
micromamba deactivate
# ##----------------------------------------------------04-02 Checkm
micromamba activate 04_checkm
cd $checkm_out || exit 1

#Run checkm on the assemblies
checkm lineage_wf -r -x fasta $all_assemblies_fasta/ $checkm_out/ --file checkm_results --tab_table -t $threads

checkm tree_qa $checkm_out/
checkm qa $checkm_out/lineage.ms $checkm_out/ -t $threads -o 2 
checkm qa $checkm_out/lineage.ms $checkm_out/ -t $threads -o 2 > $checkm_out/final_table_checkm.csv

# print the status update and output into $report
cd $home_dir
echo -e "04-02 Finished assembly quality control: Checkm. The output files are:" >> $report
for i in $checkm_out/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

micromamba deactivate
# ##----------------------------------------------------04-03 Checkm2
micromamba activate 04_checkm2

#Run checkm on the assemblies
checkm2 predict -x fasta --input $all_assemblies_fasta/ --output-directory $checkm2_out --allmodels  -t $threads 

# print the status update and output into $report
cd $home_dir
echo -e "04-03 Finished assembly quality control: Checkm2. The output files are:" >> $report
for i in $checkm2_out/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

micromamba deactivate
# ##----------------------------------------------------05_Annotate------------------------------------------- 
# ##----------------------------------------------------05_01 Prokka
micromamba activate 05_prokka

# Initialize job counter
counter=0

# Loop over all the directories in the input_dir
cd $all_assemblies_fasta || exit 1
for i in $all_assemblies_fasta/*_assembly.fasta
do  
    # Extract the base name from the file path and remove the _assembly.fasta suffix
    base=$(basename $i _assembly.fasta)

    # Run Prokka on the file
    prokka $i --outdir $prokka_out/$base --prefix $base --kingdom Bacteria --addgenes --cpus $threads_per_job_long &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

cd $prokka_out || exit 1
# Copy all gff files into one dir
find ./ -type f -name "*.gff" -exec cp {} $all_annotated_gff \;
# Copy all fsa files into one dir
find ./ -type f -name "*.fsa" -exec cp {} $all_annotated_fsa \;


# ## --------- Create a summary file containing all the information
cd $prokka_out || exit 1
# Define the output file
output_file="all_information.tsv"

# Print the header to the output file
echo -e "Tool\tSample\tcontigs\tbases\tCDS\tgene\trRNA\ttRNA\ttmRNA" > "$output_file"

# Loop through all .txt files in all subdirectories
for file in $(find . -type f -name "*.txt"); do
    sample=$(basename "$file" .txt)
    contigs=$(grep "contigs:" "$file" | cut -d' ' -f2)
    bases=$(grep "bases:" "$file" | cut -d' ' -f2)
    cds=$(grep "CDS:" "$file" | cut -d' ' -f2)
    gene=$(grep "gene:" "$file" | cut -d' ' -f2)
    rrna=$(grep "rRNA:" "$file" | cut -d' ' -f2)
    trna=$(grep "tRNA:" "$file" | cut -d' ' -f2)
    tmrna=$(grep "tmRNA:" "$file" | cut -d' ' -f2)

    # Print the values to the output file
    echo -e "unicycler\t$sample\t$contigs\t$bases\t$cds\t$gene\t$rrna\t$trna\t$tmrna" >> "$output_file"
done

# print the status update and output into $report
cd $home_dir
echo -e "05_01 Finished genome annotation: Prokka. The output files are:" >> $report
for i in $prokka_out/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

micromamba deactivate

# ##----------------------------------------------------05_02 Bakta
micromamba activate 05_bakta

for i in $all_assemblies_fasta/*.fasta; do

    # Get the base name of the file (without the path and extension)
    base=$(basename "$i" _assembly.fasta)

    # run Bakta on all samples
    bakta --db $bakta_db_dir/db --output $bakta_out/$base --prefix $base --threads $threads_per_job_long --force $i &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
        # recieve out of memory error, therefore used 8 max processes instead of $num_files_long
    if [[ $((counter%8)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

cd $bakta_out || exit 1
# Copy all gff files into one dir
find ./ -type f -name "*.gff3" -exec cp {} $bakta_gff \;
# Copy all fsa files into one dir
find ./ -type f -name "*.fna" -exec cp {} $bakta_fsa \;

# ## --------- Create a summary file containing all the information
cd $bakta_out

# Define the output file
output_file="all_information.tsv"
# Print the header to the output file
echo -e "Sample\tLength\tContig_count\tGC\tN50\tN_ratio\tcoding_density\ttRNAs\ttmRNAs\trRNAs\tncRNAs\tncRNA_regions\tCRISPR_arrays\tCDSs\tpseudogenes\thypotheticals\tsignal_peptides\tsORFs\tgaps\toriCs\toriVs\toriTs" > "$output_file"

# Loop through all .txt files in all subdirectories
for file in $(find . -type f -name "*.txt"); do
    sample=$(basename "$file" .txt)
    length=$(grep "Length:" "$file" | cut -d' ' -f2)
    contig_count=$(grep "Count:" "$file" | cut -d' ' -f2)
    gc=$(grep "GC:" "$file" | cut -d' ' -f2)
    n50=$(grep "N50:" "$file" | cut -d' ' -f2)
    n_ratio=$(grep "N ratio:" "$file" | cut -d' ' -f3)
    coding_density=$(grep "coding density:" "$file" | cut -d' ' -f3)
    trnas=$(grep "tRNAs:" "$file" | cut -d' ' -f2)
    tmrnas=$(grep "tmRNAs:" "$file" | cut -d' ' -f2)
    rrnas=$(grep "rRNAs:" "$file" | cut -d' ' -f2)
    ncrnas=$(grep "ncRNAs:" "$file" | cut -d' ' -f2)
    ncrna_regions=$(grep "ncRNA regions:" "$file" | cut -d' ' -f3)
    crispr_arrays=$(grep "CRISPR arrays:" "$file" | cut -d' ' -f3)
    cds=$(grep "CDSs:" "$file" | cut -d' ' -f2)
    pseudogenes=$(grep "pseudogenes:" "$file" | cut -d' ' -f2)
    hypotheticals=$(grep "hypotheticals:" "$file" | cut -d' ' -f2)
    signal_peptides=$(grep "signal peptides:" "$file" | cut -d' ' -f3)
    sorfs=$(grep "sORFs:" "$file" | cut -d' ' -f2)
    gaps=$(grep "gaps:" "$file" | cut -d' ' -f2)
    orics=$(grep "oriCs:" "$file" | cut -d' ' -f2)
    orivs=$(grep "oriVs:" "$file" | cut -d' ' -f2)
    orits=$(grep "oriTs:" "$file" | cut -d' ' -f2)

    # Print the values to the output file
    echo -e "$sample\t$length\t$contig_count\t$gc\t$n50\t$n_ratio\t$coding_density\t$trnas\t$tmrnas\t$rrnas\t$ncrnas\t$ncrna_regions\t$crispr_arrays\t$cds\t$pseudogenes\t$hypotheticals\t$signal_peptides\t$sorfs\t$gaps\t$orics\t$orivs\t$orits" >> "$output_file"
done

# print the status update and output into $report
cd $home_dir
echo -e "05_02 Finished genome annotation: Bakta. The output files are:" >> $report
for i in $bakta_out/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

micromamba deactivate

# ##----------------------------------------------------06_GTDB-toolkit--------------------------------------------
micromamba activate 06_gtdb

gtdbtk classify_wf --extension fasta --genome_dir $all_assemblies_fasta --out_dir $gtdb_out --cpus $threads --pplacer_cpus 1
  # pplacer issues can arise with insufficient memory, therefore reduced threads for pplacer

# print the status update and output into $report
cd $home_dir
echo -e "06 Finished classification: GTDB-tk. The output files are:" >> $report
for i in $gtdb_out/*; do
    printf -- "- %s\n" "$(basename "$i")" >> $report
done
printf "\n" >> $report

micromamba deactivate

# ##---------------------------------------------------07_Plassembler---------------------------------------------------------
micromamba activate 07_plassembler
# ##--------------------07 assembling the plasmids
cd $plassembler_out || exit 1

# Initialize job counter
counter=0

# Loop over all the long read input files, to generate base names for all samples
for file in $filtlong_out/*_filtlong.fastq.gz
do
    # Get the base name of the file
    base=$(basename "$file" _filtlong.fastq.gz)

    # Define the paired-end files
    file1=$fastp_out/${base}_1.fastq.gz
    file2=$fastp_out/${base}_2.fastq.gz
    # Define the corresponding long read (needs fastq as input, so cant use canu)
    long_read=$filtlong_out/${base}_filtlong.fastq.gz

    # Run plassembler
    plassembler run -d $plassembler_db_dir -l $long_read -o $plassembler_out/$base -1 $file1 -2 $file2 -p $base -c 500000 -t $threads_per_job_long &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# print the status update and output into $report
cd $home_dir
echo -e "07 plasmid hybrid assembly: Plassembler. The output files are:" >> $report
for i in $plassembler_out/*; do
    printf -- "- %s\n" "$i" >> $report
done
printf "\n" >> $report

micromamba deactivate
# ##-----------------------------jobinfo + deactivate mamba and unload all packages-------------------
micromamba deactivate
module purge
jobinfo 