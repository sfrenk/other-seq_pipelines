#!/usr/bin/env bash

# Maps fastq files, converters sam files to sorted + indexed bam


###############################################################################
# ARGUMENTS
###############################################################################

# Defaults
dir="."
suffix=".fastq.gz"
paired=false

usage="
	USAGE:
		This script maps fastq files and converts sam files to bam. Make sure the following modules are loaded:
		bwa
		samtools

		The script should be run from the output directory. The script should also be run with 4 processors.

		-d/--directory: directory containing fastq files. This directory should containing a subdirectory for each strain (current default: /nas02/home/s/f/sfrenk/seq/external/mmp).

		-s/--suffix: suffix of fastq files (defualt: .fastq.gz)

		-p/--paired: paired-end reads (default: false)
"

if [ -z "$1" ]; then
    echo "$usage"
    exit
fi

while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        -d|--directory)
		dir="$2"
		shift
		;;
		-s|--suffix)
		suffix="$2"
		shift
		;;
		-p|--paired)
		paired=true
		;;
    esac
shift
done


# Check directory
if [ ! -d $dir ]; then
	echo "ERROR: Invalid fastq directory path"
	exit 1
fi

# Remove trailing "/" from fastq directory if present
if [[ ${dir:(-1)} == "/" ]]; then
    dir=${dir::${#dir}-1}
fi

###############################################################################
# VARIABLES
###############################################################################

# BWA index
index="/nas02/home/s/f/sfrenk/seq/WS251/genome/bwa/genome"

#  Set the command below to display software versions used during the run
modules=$(/nas02/apps/Modules/bin/modulecmd tcsh list 2>&1)

###############################################################################
# MODULE TEST
###############################################################################

req_modules=("bwa" "samtools")

for i in ${req_modules[@]}; do
	if [[ $modules != *${i}* ]]; then
		echo "ERROR: Please load ${i}"
		exit 1
	fi
done

###############################################################################

echo "$modules"

if [ ! -d "bwa_out" ]; then
    mkdir bwa_out
fi

for file in ${dir}/*${suffix}; do
	Fbase=$(basename $file ${suffix})
    
    skipfile=false

    if [[ $paired = true ]]; then
            
        # paired end

        if [[ ${Fbase:(-2)} == "_1" ]]; then

        	base=${Fbase%_1}
        	echo "processing ${base} as paired end"
			bwa mem -t 4 $index ${dir}/${base}_1${suffix} ${dir}/${base}_2${suffix} > bwa_out/${base}.sam
		else

            # Avoid double mapping by skipping the r2 read file
                
            skipfile=true
        fi
	else
		# Single end

		base=${Fbase}
		echo "processing ${base} as single end"

		bwa mem -t 4 $index ${dir}/${base}_1.fastq.gz > bwa_out/${base}.sam
    fi

    if [[ $skipfile = false ]]; then
	# Convert to bam and remove unmapped reads (-F 4 flag)
		samtools view -bh -F 4 bwa_out/${base}.sam > bwa_out/${base}.bam
		samtools sort -o bwa_out/${base}_sorted.bam bwa_out/${base}.bam
		rm bwa_out/${base}.sam
		rm bwa_out/${base}.bam
		mv bwa_out/${base}_sorted.bam bwa_out/${base}.bam
		samtools index bwa_out/${base}.bam
	fi
done
