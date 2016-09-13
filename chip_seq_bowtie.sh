#!/usr/bin/env bash

###############################################################################
# Basic chip-seq enrichment analysis pipeline
###############################################################################

# This script produces the following output files:
# indexed, sorted bam file for each sample
# Count table for all samples combined


# Before running the script, make sure modules are loaded:
#       trim_galore
#       bowtie
#       samtools

# Check the location of the following files may have to be modified in this script:
#       bowtie index files

#       Set the command to display software versions used during the run
modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)


# NOTE: Read processing is NOT covered in this pipeline and must be performed beforehand.

###############################################################################
###############################################################################

usage="
    USAGE
       step1:   load the following modules: trim_galore bowtie samtools
       step2:   bash chip_seq_bowtie.sh [options]  

    ARGUMENTS
        -d/--dir
        Directory containing fastq files (these files should be gzipped).
        
        -p/--paired
        Use this option if fastq files contain paired-end reads. NOTE: if paired, each pair must consist of two files with the basename ending in '_1' or '_r2' depending on respective orientation.

        -t/--trim
        Trim the first n bases from the 5' end of the read (default = 0)
    "
# Set default reference

paired=false
trim=0

# Parse command line parameters

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -d|--dir)
        dir="$2"
        shift
        ;;
        -p|--paired)
        paired=true
        ;;
        -t|--trim)
        trim="$2"
        shift
        ;;
    esac
shift
done

# Remove trailing "/" from fastq directory if present

if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi


# Select bowtie index based

index="/nas02/home/s/f/sfrenk/proj/seq/WS251/genome/bowtie/genome"

# Print out loaded modules to keep a record of which software versions were used in this run

echo "$modules"

# module test

req_modules=("trim_galore" "bowtie" "samtools")

for i in ${req_modules[@]}; do
    if [[ $modules != *${i}* ]]; then
        echo "ERROR: Please load ${i}"
        exit 1
    fi
done

###############################################################################
###############################################################################

# Make directories

if [ ! -d "trimmed" ]; then
    mkdir trimmed
fi

if [ ! -d "bowtie_out" ]; then
    mkdir bowtie_out
fi

if [ ! -d "bam" ]; then
    mkdir bam
fi

# Start pipeline

echo $(date +"%m-%d-%Y_%H:%M")" Starting pipeline..."

for file in ${dir}/*.fastq.gz; do

    skipfile=false
        
    if [[ $paired = true ]]; then
            
        # paired end

        if [[ ${file:(-11)} == "_1.fastq.gz" ]]; then

            # Process and map r1 and r2 reads simultaniously

            Fbase=$(basename $file .fastq.gz)
            base=${Fbase%_1}
            
            echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} with trim_galore..."

            trim_galore --clip_R1 $trim --clip_R2 $trim --dont_gzip -o ./trimmed --paired ${dir}/${base}_1.fastq.gz ${dir}/${base}_2.fastq.gz

            echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} with bowtie..."

            # Map reads using Bowtie
            bowtie -p 4 -S -m 1 --best --strata ${index} -1 ./trimmed/${base}_1_val_1.fq -2 ./trimmed/${base}_2_val_2.fq ./bowtie_out/${base}.sam

        else

            # Avoid double mapping by skipping the r2 read file and proceding to the next step 
                
            skipfile=true
        fi
    else
            
        # Single end

        base=$(basename $file .fastq.gz)

        echo $(date +"%m-%d-%Y_%H:%M")" Trimming ${base} with trim_galore..."

        trim_galore --clip_R1 $trim --dont_gzip -o ./trimmed ${dir}/${base}.fastq.gz

        echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} with bowtie..."

        # Map reads using Bowtie
        bowtie -p 4 -S -m 1 --best --strata ${index} ./trimmed/${base}_trimmed.fq ./bowtie_out/${base}.sam
    fi     
        
    if [[ $skipfile = false ]]; then
        echo $(date +"%m-%d-%Y_%H:%M")" Mapped ${base}"

        # Convert to bam then sort

        echo $(date +"%m-%d-%Y_%H:%M")" Converting and sorting ${base}..."

        samtools view -bh -F 4 ./bowtie_out/${base}.sam > ./bam/${base}.bam

        samtools sort -o ./bam/${base}_sorted.bam ./bam/${base}.bam 

        # Need to index the sorted bam files for visualization

        echo $(date +"%m-%d-%Y_%H:%M")" indexing ${base}..."

        samtools index ./bam/${base}_sorted.bam

        rm ./bowtie_out/${base}.sam
        rm ./bam/${base}.bam
    fi

done
