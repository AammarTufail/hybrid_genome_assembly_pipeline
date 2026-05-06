#!/bin/bash

# -----------------------------------------------------------------------------
# ------------------------------ Info -----------------------------------------
# -----------------------------------------------------------------------------
# This script will install all tools and databases needed for the assembly pipeline to run

# The database downloads can be memory, cpu and time intensive so make sure you allocate as much as possible
    # tested with 28 cpus, 128 gb RAM, 48 h runtime

# Please make sure Micromamba is installed and loaded
    # https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html
        # Micromamba 1.4.2 was used for the creation of this pipeline
# Make sure Conda is also installed (only needed for the installation script)
    # GTDBtk database download did not work with micromamba

# Please change: 
    # home_dir (line 62)
    # db_home_dir (line 63)
    # Gtdb database download directory (line 200)
        # Needs to be set manually, does not accept variables as path
    # Quast: Download of the database is not working, since busco removed the old databases and quast was not yet updated to work on the new (06.05.2026)

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
# error handling
set -euo pipefail
trap 'echo "ERROR on line $LINENO"; exit 1' ERR


# ##-----------------------------Folder structure------------------------------
# Home directory
home_dir="/media/hendrik/extreme_pro1/b_assembly_pipeline_test"
db_home_dir="/media/hendrik/extreme_pro1/b_assembly_pipeline_test/db"
      #also changed the gtdb database download directory in line: 178
mkdir -p $db_home_dir
cd $home_dir || exit 1

# Databases
checkm_db_dir="$db_home_dir/checkm_db"
checkm2_db_dir="$db_home_dir/checkm2_db"
quast_db_dir="$db_home_dir/quast_db"
    # not downloadable anymore until Quast recieves a new update for the Busco database
bakta_db_dir="$db_home_dir/bakta_db"
gtdb_db_dir="$db_home_dir/gtdb_db" 
      #if path is changed, you need to update the path in line 177
plassembler_db_dir="$db_home_dir/plassembler_db"

mkdir -p $checkm_db_dir $checkm2_db_dir $quast_db_dir $bakta_db_dir $gtdb_db_dir $plassembler_db_dir


# ##-------------------------Installations--------------------------

# initialise micromamba 
eval "$(micromamba shell hook --shell bash)"
micromamba config set channel_priority strict
export MAMBA_ROOT_PREFIX=/media/hendrik/extreme_pro1/micromamba

#01- short reads qc
micromamba create -n 01_short_reads_qc python=3.11.7 -y
micromamba activate 01_short_reads_qc
micromamba install -c conda-forge -c bioconda fastqc=0.12.1 multiqc=1.19 fastp=0.23.4 -y
micromamba deactivate

#02- long reads qc
micromamba create -n 02_long_reads_qc_test python=3.10.13  -y
micromamba activate 02_long_reads_qc_test
micromamba install -c conda-forge -c bioconda longqc=1.2.0c nanofilt=2.8.0 canu=2.2 \
        filtlong=0.2.1 nanoplot=1.42.0 nanocomp=1.23.1 seqtk=1.4 -y   
pip install kaleido==0.2.1
pip install plotly==5.15.0
micromamba deactivate

#03- assembly
micromamba create -n 03_assembly python=3.10.13 -y
micromamba activate 03_assembly
micromamba install -c conda-forge -c bioconda unicycler=0.5.0  -y
micromamba install -c conda-forge gcc=13.2.0 clang=17.0.6 setuptools=69.0.3 -y
micromamba install -c conda-forge -c bioconda spades=3.15.5 racon=1.5.0 blast=2.15.0 \
        bandage=0.8.1 gapless=0.4 minimap2=2.26 seqtk=1.4 -y
micromamba install -c anaconda pandas=1.5.3 -y
micromamba deactivate

#04_assembly_qc
#04- quast
micromamba create -n 04_quast python=3.9.18 -y
micromamba activate 04_quast
micromamba install -c conda-forge -c bioconda quast=5.2.0 -y
micromamba deactivate
#04- checkm 
micromamba create -n 04_checkm python=3.9.18 -y
micromamba activate 04_checkm
micromamba install -c conda-forge -c bioconda numpy=1.26.3 matplotlib=3.8.2 pysam=0.22.0 \
         hmmer=3.4 prodigal=2.6.3 pplacer=1.1.alpha19 -y
pip3 install checkm-genome==1.2.2
micromamba deactivate

#04- checkm2
micromamba create -n 04_checkm2 python=3.12 -y
micromamba activate 04_checkm2
micromamba install -c bioconda -c conda-forge checkm2=1.1.0 --channel-priority flexible -y
micromamba deactivate

#05- Annotation
#05- Prokka
micromamba create -n 05_prokka python=3.12.1 -y
micromamba activate 05_prokka
micromamba install -c conda-forge -c bioconda prokka=1.14.6 -y
micromamba deactivate

#05- Bakta
micromamba create -n 05_bakta \
    -c conda-forge \
    -c bioconda \
    python=3.9 \
    bakta=1.9.3 \
    pyhmmer=0.10.15 \
    -y

#06- GTDB-Tk
micromamba create -n 06_gtdb python=3.8.15 -y
micromamba activate 06_gtdb
micromamba install -c conda-forge -c bioconda gtdbtk=2.1.1 numpy=1.23.1 -y
micromamba deactivate

#07- Plassembler
micromamba create -n 07_plassembler_new python=3.9.18 -y
micromamba activate 07_plassembler_new
micromamba install -c conda-forge -c bioconda plassembler=1.8.2 -y
micromamba deactivate

# # ##-------------------------Download Databases--------------------------
#Checkm
micromamba activate 04_checkm
cd $checkm_db_dir
#download the database file
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
#extract the database file
tar xvzf checkm_data_2015_01_16.tar.gz 
#create the database
export CHECKM_DATA_PATH=$checkm_db_dir
checkm data setRoot $checkm_db_dir
micromamba deactivate

#Checkm2
micromamba activate 04_checkm2
checkm2 database --download --path $checkm2_db_dir 
micromamba deactivate

# #Quast
# micromamba activate 04_quast
# #Download of the database is not working, since busco removed the old databases and quast was not updated to work on the new db yet (04.02.2024)
# cd $HOME/.micromamba/envs/04_quast/lib/python3.9/site-packages/quast_libs/augustus3.2.3
# wget http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.2.3.tar.gz
# tar xvzf augustus-3.2.3.tar.gz
# micromamba deactivate

#Bakta
micromamba activate 05_bakta
bakta_db download --output $bakta_db_dir --type full
micromamba deactivate

#Gtdb-tk
 micromamba activate 06_gtdb 
 cd $gtdb_db_dir || exit 1
# Download the database file
 wget -c --tries=5 https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
ls gtdbtk_r220_data.tar.gz || { echo "Error: gtdbtk_r220_data.tar.gz not found. Exit script"; exit 1; }
tar -xvzf gtdbtk_r220_data.tar.gz
#set the db path
micromamba env config vars set GTDBTK_DATA_PATH="/media/hendrik/extreme_pro1/b_assembly_pipeline_test/db/gtdb_db/release220"
export GTDBTK_DATA_PATH="$gtdb_db_dir/release220"
micromamba deactivate

#Plassembler
micromamba activate 07_plassembler
plassembler download -d $plassembler_db_dir
micromamba deactivate
