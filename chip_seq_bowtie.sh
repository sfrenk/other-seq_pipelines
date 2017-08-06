#!/usr/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 1-0

module load bowtie samtools deeptools
module list

###############################################################################
# Basic chip-seq enrichment analysis pipeline
###############################################################################

# This script produces the following output files for each sample:
#   indexed, sorted bam file
#   normalized bedgraph and BigWig files


# NOTE: Read trimming/quality filtering is NOT covered in this pipeline (with the exception of 5' and 3' hard trimming) and must be performed beforehand using trim.sh.

###############################################################################
###############################################################################

usage="
    USAGE
      sbatch chip_seq_bowtie.sh [options]  

    ARGUMENTS
        -d/--dir
        Directory containing fastq files (these files should be gzipped).
        
        -p/--paired
        Use this option if fastq files contain paired-end reads. NOTE: if paired, each pair must consist of two files with the basename ending in '_1' or '_2' depending on respective orientation.

        -m/--multi
        How to deal with reads with more than one best alignment:

            random (default: randomly assign the read to one of the locations
            discard: discard the read
            all: report all alignments

        -l/--left
        Bases to trim from left (5') end

        -r/--right
        Bases to trim from right (3') end

        -e/--extend
        Extend reads (if using single end reads, supply the mean fragment length. For paired end, deeptools calculates this automatically)
    "
# Set default reference

paired=false
multi_option="random"
left=0
right=0
frag_flag=""
sam_flag=""

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
        -m|--multi)
        multi_option="$2"
        shift
        ;;
        -l|--left)
        left="$2"
        shift
        ;;
        -r|--right)
        right="$2"
        shift
        ;;
        -e|--extend)
        if [[ $2 =~ '^[0-9]+$' ]] ; then
            frag_flag="--extendReads $2"
        else:
            frag_flag="--extendReads"
        fi
        shift
        ;;
    esac
    shift
done


# bowtie index

index="/nas02/home/s/f/sfrenk/proj/seq/WS251/genome/bowtie/genome"

# parse multi_option

case $multi_option in
    "discard")
    multi_flag="-m 1"
    ;;
    "random")
    multi_flag="-M 1"
    ;;
    "all")
    multi_flag="-a"
    ;;
esac

# Set the --samFlag option for bamCoverage so that only foward reads are counted

if [[ $paired = true ]]; then
    sam_flag="--samFlagInclude 64"
fi

# Print run parameters to file

if [[ ! -f run_parameters.txt ]]; then
    printf "$(date +"%m-%d-%Y_%H:%M")\n\nPARAMETERS\n\nPipeline: chip_seq_bowtie\n\nsample directory: ${dir}\nmulti: ${multi_flag}\npaired end: ${paired}\nleft_trim: ${left}\nright_trim: ${right}\nread_extension: ${frag_flag}\n\n" > run_parameters.txt

    module list &>> run_parameters.txt
    printf "\n" >> run_parameters.txt
fi

###############################################################################
###############################################################################

# Make directories

if [ ! -d "bowtie_out" ]; then
    mkdir bowtie_out
fi

if [ ! -d "bam" ]; then
    mkdir bam
fi

if [ ! -d "cov" ]; then
    mkdir cov
fi

# Start pipeline

echo $(date +"%m-%d-%Y_%H:%M")" Starting pipeline..."

shopt -s nullglob

files=(${dir}/*.fastq.gz)

for file in ${files[@]}; do

    skipfile=false
        
    if [[ $paired = true ]]; then
            
        # paired end

        if [[ ${file:(-11)} == "_1.fastq.gz" ]]; then

            # Process and map r1 and r2 reads simultaniously

            Fbase=$(basename $file .fastq.gz)
            base=${Fbase%_1}
            
            echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} with bowtie..."

            # Map reads using Bowtie

            gunzip -c ${dir}/${base}_1.fastq.gz > ${base}_1.fastq
            gunzip -c ${dir}/${base}_2.fastq.gz > ${base}_2.fastq 

            bowtie -p $SLURM_NTASKS -S $multi_flag --chunkmbs 200 --best --strata -5 $left -3 $right ${index} -1 ${base}_1.fastq -2 ${base}_2.fastq ./bowtie_out/${base}.sam

            rm ${base}_*.fastq

        else

            # Avoid double mapping by skipping the r2 read file and proceding to the next step 
                
            skipfile=true
        fi
    else
            
        # Single end

        base=$(basename $file .fastq.gz)

        echo $(date +"%m-%d-%Y_%H:%M")" Mapping ${base} with bowtie..."

        # Map reads using Bowtie
        gunzip -c $file > ${base}.fastq

        bowtie -p $SLURM_NTASKS -S $multi_flag --chunkmbs 200 --best --strata -5 $left -3 $right ${index} ${base}.fastq ./bowtie_out/${base}.sam

        rm ${base}.fastq
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

        total_mapped="$(samtools view -c ./bam/${base}_sorted.bam)"
        printf "${base}\t${total_mapped}\n" >> total_mapped_reads.txt

        # Make normalized genome track files
        printf "Making normalized bg and bw files...\n"

        # For visualization
        bamCoverage -b ./bam/${base}_sorted.bam -o ./cov/${base}.bw -bs 10 -p $SLURM_NTASKS --normalizeUsingRPKM $frag_flag $sam_flag

        # For quantification
        bamCoverage -b ./bam/${base}_sorted.bam -o ./cov/${base}.1kb.bg -of bedgraph -bs 200 -p $SLURM_NTASKS --normalizeUsingRPKM $frag_flag $sam_flag
    fi

done
