#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=23:55:00
#SBATCH --job-name=tool_test
#SBATCH --output=tool_test.out
#SBATCH --error=tool_test.err
#SBATCH --partition=base

# -----------------------------------------------------------------------------
# ------------------------------ Info -----------------------------------------
# -----------------------------------------------------------------------------
# Installation, commands and quality analysis of every genome assembly tool we tested 
    # Please change: 
        # home_dir (line 78)
        # Micromamba activation (line 67-69)
        # proxy settings for HPC without direct internet connection (line 71-74)
# Already processed reads are required this script 
  # workflow found in assembly_pipeline.sh
  # adjust read_dir path to them (line 92-99)
# Installation and database download for all assembly_qc (Quast, Checkm, Checkm2, Prokka, Bakta) tools is required
  # found in installation_pipeline.sh
  # Quality assessmant tools start at line 790


#--------------Tools:
## Assembly
# SKESA (https://github.com/ncbi/SKESA)
    # Short read assembly
# SPAdes (https://github.com/ablab/spades)
    # Short read assembly
# Wtdbg2 (https://github.com/ruanjue/wtdbg2)
    # Long read assembly
# Miniasm (https://github.com/lh3/miniasm)
    # Long read assembly
# Flye (https://github.com/lh3/miniasm)
    # Long read assembly
# Hybridspades  (https://github.com/ablab/spades ; https://ablab.github.io/spades/hybrid.html)
  # Hybrdid assembly
# Unicycler (https://github.com/rrwick/Unicycler
    # Hybrid genome assembly from short and long reads
# Canu (https://github.com/marbl/canu)
    # Error correction of ONT long reads
# Racon (https://ablab.github.io/spades/hybrid.html)
  # Asembly polishing/correcting with short and long reads

## Quality control
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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Load modules
module load micromamba/1.4.2
export MAMBA_ROOT_PREFIX=$HOME/.micromamba
eval "$(micromamba shell hook --shell=bash)"

#set proxy environment for HPC
export http_proxy=http://relay:3128
export https_proxy=http://relay:3128
export ftp_proxy=http://relay:3128

# ##--------------------------------Set directory paths-----------------------------
# Home directory
home_dir="/path/to/your/assembly_pipeline_directory"
mkdir -p $home_dir
cd $home_dir || exit 1

# Databases
db_home_dir="/path/to/your/databases_directory"
checkm_db_dir="$db_home_dir/checkm_db"
checkm2_db_dir="$db_home_dir/checkm2_db"
# quast_db_dir="$db_home_dir/quast_db"
    # not downloadable anymore until Quast recieves a new update for the Busco database
bakta_db_dir="$db_home_dir/bakta_db"

# Read input    coot
read_dir="/path/to/your/reads" 
  # 01-Trimmed short reads
  fastp_out="$read_dir/01_short_reads"
  # 02-Trimmed long reads
  filtlong_out_fasta="$read_dir/02_long_reads"
  # 03-Canu corrected long reads
    # only needed for unicycler and racon
  canu_corrected_reads="$read_dir/03_corrected_long_reads"

#---------------------------------------------------------------------------------------
# ##-----Set the amount of Ram and CPUs available per task running parallel-------------
#---------------------------------------------------------------------------------------

# define the maximum amount of threads available
threads=$(nproc)

# Total number of input files
num_files_short=$(ls $fastp_out/* | wc -l)
num_files_long=$(ls $filtlong_out_fasta/* | wc -l)

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

#---------------------------------------------------------------------------------------
# ##-----------------------------------Short read assemblies----------------------------
#---------------------------------------------------------------------------------------

#----------------------------- Skesa
# ##Installation
micromamba create -n test_skesa -y
micromamba activate test_skesa
micromamba install bioconda::skesa=2.5.1 -y

# ##Running
skesa_out="$home_dir/test_skesa"
mkdir -p $skesa_out
cd $skesa_out 

# Initialize job counter
counter=0
# Loop over all cleaned short reads
for file1 in $fastp_out/*_1.fastq.gz
do  
    # Extract basename of the file
    base=$(basename $file1 _1.fastq.gz)
      
    # Define paired file
    file2=${file1/_1.fastq.gz/_2.fastq.gz}
    
    #create a output folder for the sample
    mkdir -p $base

    # Run wgtdb2 on the file
    skesa --reads $file1,$file2 --cores $threads_per_job_short > $base/${base}.skesa.fa &
    counter=$((counter+1))
    
    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_short)) -eq 0 ]]; then
        wait
        counter=0
    fi   
done
# Wait for all background jobs to complete
wait

# Copy all the final assemblies into one folder
all_assemblies_skesa="$skesa_out/all_assemblies"
mkdir -p $all_assemblies_skesa
# Loop over all the directories in the input_dir
for dir in $skesa_out/*/
do    
    # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_skesa" ]; then
  
      # Extract the subdirectory name (without the trailing slash)
      base=$(basename "$dir")

      # Copy and rename the file to the target directory
      cp $dir/*.skesa.fa $all_assemblies_skesa/${base}.fasta
  fi
done


micromamba deactivate

#----------------------------- Spades
# ##Installation
micromamba create -n test_spades -y
micromamba activate test_spades
micromamba install -c bioconda spades=3.15.5 -y

# ##Running
spades_out="$home_dir/test_spades"
mkdir -p $spades_out
cd $spades_out

# Initialize job counter
counter=0
# Loop over all cleaned short reads
for i in $fastp_out/*_1.fastq.gz
do
    
    # Extract basename of the file
    base=$(basename "$i" _1.fastq.gz)

    # Define the paired-end files
    file1=$fastp_out/${base}_1.fastq.gz
    file2=$fastp_out/${base}_2.fastq.gz


    # Run Unicycler on the paired-end files and the long reads
    spades.py -1 $file1 -2 $file2 -o $base --isolate -t $threads_per_job_long &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_short)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# Copy all the final assemblies into one folder
all_assemblies_spades="$spades_out/all_assemblies"
mkdir -p $all_assemblies_spades
# Loop over all the directories in the input_dir
for dir in $spades_out/*
do    
    # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_spades" ]; then
  
      # Extract the subdirectory name (without the trailing slash)
      base=$(basename "$dir")
      
      # Copy and rename the file to the target directory
      cp $dir/scaffolds.fasta $all_assemblies_spades/${base}.fasta
  fi
done

micromamba deactivate

#---------------------------------------------------------------------------------------
# ##-----------------------------------Long read assemblies-----------------------------
#---------------------------------------------------------------------------------------

#------------------------------ Wtdbg2
# ##Installation
micromamba create -n test_wtdbg -y
micromamba activate test_wtdbg
micromamba install bioconda::wtdbg=2.5 -y

# ##Running
wtdbg2_out="$home_dir/test_wtdbg2"
mkdir -p $wtdbg2_out
cd $wtdbg2_out

# #Contig assembly
# Initialize job counter
counter=0
# Loop over all cleaned long reads
for i in $filtlong_out_fasta/*fasta.gz
do  
    # Extract basename of the file
    base=$(basename $i _filtlong.fasta.gz)
    
    #create a output folder for the sample
    mkdir -p $base

    # Run wgtdb2 on the file
    wtdbg2 -x ont -g 5m -t $threads_per_job_long -i $i -fo $base/$base &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# #Create the consensus of all contigs
cd $wtdbg2_out
# Initialize job counter
counter=0
# Loop over all the directories in the input_dir
for i in */
do  
    cd $i
    # Extract the base name (remove trailing slash)
    base=$(basename $i /)
    
    # Run wgtdb2 on the file
    wtpoa-cns -t $threads_per_job_long -i ${base}.ctg.lay.gz -fo ${base}.ctg.fa &
    counter=$((counter+1))
    
    cd ../

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# #Copy all the final assemblies into one folder
all_assemblies_wtdbg2="$wtdbg2_out/all_assemblies"
mkdir -p $all_assemblies_wtdbg2
# Loop over all the directories in the input_dir
for dir in $wtdbg2_out/*
do    
    # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_wtdbg2" ]; then
  
      # Extract the subdirectory name (without the trailing slash)
      base=$(basename "$dir")

      # Copy and rename the file to the target directory
      cp $dir/*.ctg.fa $all_assemblies_wtdbg2/${base}.fasta
  fi
done


micromamba deactivate

#------------------------------ Miniasm
# ##Installation
micromamba create -n test_miniasm -y
micromamba activate test_miniasm
micromamba install -c bioconda miniasm=0.3 minimap2=2.28 gfatools=0.5 -y

# ##Running
miniasm_out="$home_dir/test_miniasm"
mkdir -p $miniasm_out
cd $miniasm_out

# #Read mapping
# Initialize job counter
counter=0
# Loop over all cleaned long reads
for i in $filtlong_out_fasta/*fasta.gz
do  
    # Extract basename of the file
    base=$(basename $i _filtlong.fasta.gz)
    
    #create a output folder for the sample
    mkdir -p $base

    # Run wgtdb2 on the file
    minimap2 -x ava-ont -t $threads_per_job_long $i $i | gzip -1 > $base/${base}.paf.gz &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# Create the consensus of all contigs------------------------------
cd $miniasm_out
# Initialize job counter
counter=0
# Loop over all the directories in the input_dir
for dir in $miniasm_out/*
do  
    # Extract basename of the directory
    base=$(basename $dir)
    
    cd $base
    echo -e  "miniasm running in folder $base : $dir" 
    
    # Run miniasm on the file
    miniasm -f $filtlong_out_fasta/${base}* ${base}.paf.gz > ${base}.gfa 
    
    cd ../

done

# #Convert .gfa to .fasta and place all assemblies into one folder
all_assemblies_miniasm="$miniasm_out/all_assemblies"
mkdir -p $all_assemblies_miniasm
# Loop over all the directories in the input_dir
for dir in $miniasm_out/*
do  
  # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_miniasm" ]; then
  
    # Extract the subdirectory name (without the trailing slash)
    base=$(basename "$dir")
    
    # convert the gfa to fasta and move&rename to the target directory
    gfatools gfa2fa $dir/*.gfa > $all_assemblies_miniasm/${base}.fasta
  fi
done


micromamba deactivate

#------------------------------ Flye
# ##Installation
micromamba create -n test_flye -y
micromamba activate test_flye
micromamba install -c bioconda flye=2.9.4 -y

flye_out="$home_dir/test_flye"
mkdir -p $flye_out
cd $flye_out

# #assemble the reads into contigs
# Initialize job counter
counter=0
# Loop over all cleaned long reads
for i in $filtlong_out_fasta/*fasta.gz
do  
    # Extract basename of the file
    base=$(basename $i _filtlong.fasta.gz)
    
    #create a output folder for the sample
    mkdir -p $base

    # Run wgtdb2 on the file
    flye --nano-corr $i --out-dir $base -g 5m --threads $threads_per_job_long &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# #Copy all the final assemblies into one folder
all_assemblies_flye="$flye_out/all_assemblies"
mkdir -p $all_assemblies_flye
# Loop over all the directories in the input_dir
for dir in $flye_out/*
do 
  # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_$flye" ]; then
  
    # Extract the subdirectory name (without the trailing slash)
    base=$(basename "$dir")
    
    # convert the gfa to fasta and move&rename to the target directory
    cp $dir/assembly.fasta $all_assemblies_flye/${base}.fasta
  fi
done

micromamba deactivate

#---------------------------------------------------------------------------------------
# ##-----------------------------------Hybrid assemblies--------------------------------
#---------------------------------------------------------------------------------------

#------------------------------ Hybridspades
# ##Installation
micromamba create -n test_spades -y
micromamba activate test_spades
micromamba install -c bioconda spades=3.15.5 -y

# ##Running
hybridspades_out="$home_dir/test_hybridspades"
mkdir -p $hybridspades_out
cd $hybridspades_out

# Initialize job counter
counter=0
# Loop over all cleaned long reads
for i in $filtlong_out_fasta/*fasta.gz
do  
    # Extract basename of the file
    base=$(basename $i _filtlong.fasta.gz)

    # Define the paired-end files
    file1=$fastp_out/${base}_1.fastq.gz
    file2=$fastp_out/${base}_2.fastq.gz
    
    #create a output folder for the sample
    mkdir -p $base

    # Run whybridspades on the file
    spades.py --isolate -1 $file1 -2 $file2 \
              --nanopore $i \
              -o $base -t $threads_per_job_long &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# #Copy all the final assemblies into one folder
all_assemblies_hybridspades="$hybridspades_out/all_assemblies"
mkdir -p $all_assemblies_hybridspades
# Loop over all the directories in the input_dir
for dir in $hybridspades_out/*
do    
    # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_hybridspades" ]; then
  
      # Extract the subdirectory name (without the trailing slash)
      base=$(basename "$dir")

      # Copy and rename the file to the target directory
      cp $dir/scaffolds.fasta $all_assemblies_hybridspades/${base}.fasta
  fi
done

micromamba deactivate

#------------------------------ Unicycler
# ##Installation
micromamba create -n test_unicycler python=3.10.13 -y
micromamba activate test_unicycler
micromamba install -c bioconda unicycler=0.5.0  -y
micromamba install -c conda-forge gcc=13.2.0 clang=17.0.6 setuptools=69.0.3 -y
micromamba install -c bioconda spades=3.15.5 racon=1.5.0 blast=2.15.0 \
        bandage=0.8.1 gapless=0.4 minimap2=2.26 seqtk=1.4 -y
micromamba install -c anaconda pandas=1.5.3 -y

# ##Running
unicycler_out="$home_dir/test_unicycler"
mkdir -p $unicycler_out
cd $unicycler_out

# Initialize job counter
counter=0
# Loop over all cleaned long reads
for i in $filtlong_out_fasta/*fasta.gz
do  
    # Extract basename of the file
    base=$(basename $i _filtlong.fasta.gz)

    # Define the paired-end files
    file1=$fastp_out/${base}_1.fastq.gz
    file2=$fastp_out/${base}_2.fastq.gz

    # Run Unicycler on the paired-end files and the long reads
    unicycler -1 $file1 -2 $file2 -l $i -o $unicycler_out/$base -t $threads_per_job_long &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# #Copy all the final assemblies into one folder
all_assemblies_unicycler="$unicycler_out/all_assemblies"
mkdir -p $all_assemblies_unicycler
# Loop over all the directories in the input_dir
for dir in $unicycler_out/*
do 
    # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_unicycler" ]; then
  
      # Extract the subdirectory name (without the trailing slash)
      base=$(basename "$dir")

      # Copy and rename the file to the target directory
      cp $dir/assembly.fasta $all_assemblies_unicycler/${base}.fasta
  fi
done

micromamba deactivate

#------------------------------ Unicycler with Canu
# ##Installation
micromamba activate test_unicycler


# ##Running
unicycler_canu_out="$home_dir/test_unicycler_canu"
mkdir -p $unicycler_canu_out
cd $unicycler_canu_out

# Initialize job counter
counter=0
# Loop over all cleaned long reads
for i in $canu_corrected_reads/*fasta.gz
do  
    # Extract basename of the file
    base=$(basename $i .correctedReads.fasta.gz)

    # Define the paired-end files
    file1=$fastp_out/${base}_1.fastq.gz
    file2=$fastp_out/${base}_2.fastq.gz

    # Run Unicycler on the paired-end files and the long reads
    unicycler -1 $file1 -2 $file2 -l $i -o $unicycler_canu_out/$base -t $threads_per_job_long &
    counter=$((counter+1))

    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# #Copy all the final assemblies into one folder
all_assemblies_unicycler_canu="$unicycler_canu_out/all_assemblies"
mkdir -p $all_assemblies_unicycler_canu
# Loop over all the directories in the input_dir
for dir in $unicycler_canu_out/*
do  
    # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_unicycler_canu" ]; then
  
      # Extract the subdirectory name (without the trailing slash)
      base=$(basename "$dir")

      # Copy and rename the file to the target directory
      cp $dir/assembly.fasta $all_assemblies_unicycler_canu/${base}.fasta
  fi
done

micromamba deactivate

#---------------------------------------------------------------------------------------
# ##-----------------------------------Assembly polishing-------------------------------
#---------------------------------------------------------------------------------------
# ------------------------------- Racon
# polished the canu corrected unicycler assembly (best results so far)
  # VERY memory intensive 
      # did not run with less than 64G, even without any multi-threading
      
# ##Installation
micromamba create -n test_racon -y
micromamba activate test_racon
micromamba install -c bioconda racon=1.5.0 minimap2=2.28 bwa=0.7.18 samtools=1.20  -y

# ##Running
racon_out="$home_dir/test_unicycler_canu_racon"
mkdir -p $racon_out
cd $racon_out


#--------------long read polishing

#align/map the long reads to the assembly
# Initialize job counter
counter=0
# Loop over all assemblies
for i in $all_assemblies_unicycler_canu/*.fasta
do  
    # Extract the base name from the file 
    base=$(basename $i .fasta)   
    #create a output folder for the sample
    mkdir -p $base
    cd $base

    # Run the command
    minimap2 -x map-ont -t $threads_per_job_long $i $canu_corrected_reads/${base}.correctedReads.fasta.gz > ${base}_lr_mapped.paf &
    counter=$((counter+1))

    cd ../    
    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# polish the assembly with the aligned long reads
# Initialize job counter
counter=0
# Loop over all assemblies
for i in $all_assemblies_unicycler_canu/*.fasta
do  
    # Extract the base name from the file 
    base=$(basename $i .fasta)         
    #go to output folder for the sample
    cd $base

    # Run the command
    racon -t $threads_per_job_long $canu_corrected_reads/${base}.correctedReads.fasta.gz ${base}_lr_mapped.paf $i > ${base}_lr_polished.fasta &
    counter=$((counter+1))

    cd ../ 
    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait


#-------------------------------short read polishing

# index the assemblies polished with the long reads 
# Initialize job counter
counter=0
# Loop over all assemblies
for i in $all_assemblies_unicycler_canu/*.fasta
do  
    # Extract the base name from the file 
    base=$(basename $i .fasta)           
    #go to output folder for the sample
    cd $base

    # Run the command
    bwa index ${base}_lr_polished.fasta &
    counter=$((counter+1))
    
    cd ../
    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# align the short reads to the lr polished assemblies
# Initialize job counter
counter=0
# Loop over all assemblies
for i in $all_assemblies_unicycler_canu/*.fasta
do  
    # Extract the base name from the file 
    base=$(basename $i .fasta)           
    #go to output folder for the sample
    cd $base

    # Run the command
    bwa mem -t $threads_per_job_long ${base}_lr_polished.fasta $fastp_out/${base}_1.fastq.gz $fastp_out/${base}_2.fastq.gz > ${base}_sr_mapped.sam &
    counter=$((counter+1))

    cd ../     
    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# sort the sam file (can`t use .bam for racon)
# Initialize job counter
counter=0
# Loop over all assemblies
for i in $all_assemblies_unicycler_canu/*.fasta
do  
    # Extract the base name from the file 
    base=$(basename $i .fasta)      
    #go to output folder for the sample
    cd $base

    # Run the command
    samtools sort -O SAM -o ${base}_sr_sorted.sam ${base}_sr_mapped.sam &
    counter=$((counter+1))

    cd ../     
    # If the number of jobs has reached the maximum, wait for all jobs to finish before starting more
    if [[ $((counter%num_files_long)) -eq 0 ]]; then
        wait
        counter=0
    fi
done
# Wait for all background jobs to finish
wait

# polish the lr polished assembly with aligned short reads
    # Run without parallelisation, for more ram per sample
# Initialize job counter
#counter=0
# Loop over all assemblies
for i in $all_assemblies_unicycler_canu/*.fasta
do  
    # Extract the base name from the file 
    base=$(basename $i .fasta) 
    #go to output folder for the sample
    cd $base
   
    # concatenate the foreward and reverse reads into one file for racoon to use.
    cat $fastp_out/${base}_1.fastq.gz $fastp_out/${base}_2.fastq.gz > ${base}_combined.fastq.gz

    # run racon
    racon -t $threads_per_job_long -u ${base}_combined.fastq.gz ${base}_sr_sorted.sam ${base}_lr_polished.fasta > ${base}_final_polished.fasta 
    
    cd ../    
done
# Wait for all background jobs to finish
wait

# #Copy all the final assemblies into one folder
all_assemblies_racon="$racon_out/all_assemblies"
mkdir -p $all_assemblies_racon
# Loop over all the directories in the input_dir
for dir in $racon_out/*
do  
    # Skip the 'all_assemblies' directory
  if [ "$dir" != "$all_assemblies_racon" ]; then
  
      # Extract the subdirectory name (without the trailing slash)
      base=$(basename "$dir")

      # Copy and rename the file to the target directory
      cp $dir/*_final_polished.fasta $all_assemblies_racon/${base}.fasta
  fi
done

micromamba deactivate

#---------------------------------------------------------------------------------------
# ##-----------------------------------Quality assessment-------------------------------
#---------------------------------------------------------------------------------------
# All 5 tools used to assess the quality will loop over every file in each "all_assemblies" dir of each tested tool
  # results will be placed in one dir, with subdirs for each tested tool

home_dir="/work_beegfs/sunam216/hendrik/paper_script_test/tools_test"

# Loop over all subdirectories that contain an 'all_assemblies' directory
for dir in "$home_dir"/test_*; do
  
    if [ -d "$dir/all_assemblies/" ]; then
        # save the name of the tool (remove 'test_' from the directory name)
        tool=$(basename "$dir" | sed 's/^test_//')
        # Define output directories specific to each subdirectory
        quast_out="$home_dir/quality_assessment/quast/$tool"
        checkm_out="$home_dir/quality_assessment/checkm/$tool"
        checkm2_out="$home_dir/quality_assessment/checkm2/$tool"
        prokka_out="$home_dir/quality_assessment/prokka/$tool"
        bakta_out="$home_dir/quality_assessment/bakta/$tool"
        mkdir -p $quast_out $checkm_out $checkm2_out $prokka_out $bakta_out 

        #----------------------04-01 Quast
        micromamba activate 04_quast
        cd $dir/all_assemblies
        quast.py -o "$quast_out" *.fasta --circos --glimmer --rna-finding --conserved-genes-finding --report-all-metrics --use-all-alignments -t $threads
        micromamba deactivate

        # ----------------------04-02 Checkm
        micromamba activate 04_checkm
        
        cd $checkm_out
        checkm lineage_wf -r -x fasta $dir/all_assemblies/ $checkm_out/ --file checkm_results --tab_table -t $threads
        checkm tree_qa $checkm_out/
        checkm qa $checkm_out/lineage.ms $checkm_out/ -t $threads -o 2 > $checkm_out/final_table_checkm.csv
        micromamba deactivate

        # ----------------------04-03 Checkm2
        micromamba activate 04_checkm2
        checkm2 predict -x fasta --input $dir/all_assemblies/ --output-directory $checkm2_out --allmodels -t $threads
        micromamba deactivate

        # ----------------------05-01 Prokka
        micromamba activate 05_prokka
        cd $dir/all_assemblies
        counter=0
        for i in ./*.fasta; do
            base=$(basename "$i" .fasta)
            prokka $i --outdir $prokka_out/$base --prefix $base --kingdom Bacteria --addgenes --cpus $threads_per_job_long &
            counter=$((counter+1))
            if [[ $((counter%num_files_long)) -eq 0 ]]; then
                wait
                counter=0
            fi
        done
        wait
        micromamba deactivate

        # -----------------------05-02 Bakta
        micromamba activate 05_bakta
        cd $dir/all_assemblies
        counter=0
        for i in ./*.fasta; do
            base=$(basename "$i" _assembly.fasta)
            bakta --db $bakta_db_dir/db --output $bakta_out/$base --prefix $base --threads $threads_per_job_long --force $i &
            counter=$((counter+1))
            if [[ $((counter%num_files_long)) -eq 0 ]]; then
                wait
                counter=0
            fi
        done
        wait       
        micromamba deactivate
  fi
done

wait

# -------------------- Concatenate the results of each quality control tool
# Create summary file for Quast
cd $quast_out
cd ../
output_quast="all_information_quast.tsv"
for file in $(find . -type f -name "report.txt"); do
    tool=$(basename $(dirname "$file") | sed 's/^test_//')
    # Check if the output file exists. If not:
    if [ ! -f $output_quast ]; then
      # Skip the first two lines and process from the "Assembly" header
      tail -n +3 $file | while read -r line
      do
        # Append the tool name to each row and write to the output file
        echo -e "$tool\t$line" >> "$output_quast"
        echo -e "creating file with header for $tool, containing $line in file $output_quast" # testing
      done
    else 
      # Skip the first two lines and process from the "Assembly" header
      tail -n +4 $file | while read -r line
      do
        # Append the tool name to each row and write to the output file
        echo -e "$tool\t$line" >> "$output_quast"
      done
    fi
done

# Create summary file for checkm
cd $checkm_out
cd ../
output_checkm="all_information_checkm.tsv"
for file in $(find . -type f -name "final_table_checkm.csv"); do
    tool=$(basename $(dirname "$file") | sed 's/^test_//')
    # Check if the output file exists. If not:
    if [ ! -f $output_checkm ]; then
      # copy the entire table including the header
      tail -n +9 $file | while read -r line
      do
        # Append the tool name to each row and write to the output file
        echo -e "$tool\t$line" >> "$output_checkm"
      done
    else 
      # Skip the first line and process without the header
      tail -n +10 $file | while read -r line
      do
        # Append the tool name to each row and write to the output file
        echo -e "$tool\t$line" >> "$output_checkm"
      done
    fi
done

# Create summary file for checkm2
cd $checkm2_out
cd ../
output_checkm2="all_information_checkm2.tsv"
for file in $(find . -type f -name "quality_report.tsv"); do
    tool=$(basename $(dirname "$file") | sed 's/^test_//')
    # Check if the output file exists. If not:
    if [ ! -f $output_checkm2 ]; then
      # copy the entire table including the header
      tail $file | while read -r line
      do
        # Append the tool name to each row and write to the output file
        echo -e "$tool\t$line" >> "$output_checkm2"
      done
    else 
      # Skip the first line and process without the header
      tail -n +2 $file | while read -r line
      do
        # Append the tool name to each row and write to the output file
        echo -e "$tool\t$line" >> "$output_checkm2"
      done
    fi
done


# Create summary file for Prokka
cd $prokka_out
cd ../
output_prokka="all_information_prokka.tsv"
echo -e "Tool\tSample\tcontigs\tbases\tCDS\tgene\trRNA\ttRNA\ttmRNA" > "$output_prokka"
for file in $(find . -type f -name "*.txt"); do
    tool=$(basename $(dirname "$file") | sed 's/^test_//')
    sample=$(basename "$file" .txt)
    contigs=$(grep "contigs:" "$file" | cut -d' ' -f2)
    bases=$(grep "bases:" "$file" | cut -d' ' -f2)
    cds=$(grep "CDS:" "$file" | cut -d' ' -f2)
    gene=$(grep "gene:" "$file" | cut -d' ' -f2)
    rrna=$(grep "rRNA:" "$file" | cut -d' ' -f2)
    trna=$(grep "tRNA:" "$file" | cut -d' ' -f2)
    tmrna=$(grep "tmRNA:" "$file" | cut -d' ' -f2)
    echo -e "$tool\t$sample\t$contigs\t$bases\t$cds\t$gene\t$rrna\t$trna\t$tmrna" >> "$output_prokka"
done
  

# Create summary file for Bakta
cd $bakta_out
cd ../
output_bakta="all_information_bakta.tsv"
echo -e "Tool\tSample\tLength\tContig_count\tGC\tN50\tN_ratio\tcoding_density\ttRNAs\ttmRNAs\trRNAs\tncRNAs\tncRNA_regions\tCRISPR_arrays\tCDSs\tpseudogenes\thypotheticals\tsignal_peptides\tsORFs\tgaps\toriCs\toriVs\toriTs" > "$output_bakta"
for file in $(find . -type f -name "*.txt"); do
    tool=$(basename $(dirname "$file") | sed 's/^test_//')
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
    echo -e "$tool\t$sample\t$length\t$contig_count\t$gc\t$n50\t$n_ratio\t$coding_density\t$trnas\t$tmrnas\t$rrnas\t$ncrnas\t$ncrna_regions\t$crispr_arrays\t$cds\t$pseudogenes\t$hypotheticals\t$signal_peptides\t$sorfs\t$gaps\t$orics\t$orivs\t$orits" >> "$output_bakta"
done

# ##-----------------------------jobinfo + deactivate mamba and unload all packages-------------------
micromamba deactivate
module purge
jobinfo 