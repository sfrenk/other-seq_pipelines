import glob
import re
import sys

############################    DESCRIPTION    ##############################

# Call Variants using GATK
# NOTE: use setup_dir.sh -d <directory_name> -p gatk to set up the working directory

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = ""
EXTENSION = ""
PAIRED = False

# Trimming parameters
ADAPTERS="/nas/longleaf/apps/bbmap/37.62/bbmap/resources/adapters.fa"

# Mapping parameters
BWA_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/bwa/genome"

###############################################################################

SAMPLE_FILES = glob.glob(BASEDIR + "/*" + EXTENSION)
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + EXTENSION, x).group(1) for x in SAMPLE_FILES ]

# Paired end files must end in _1/_2 where 1 and 2 denote forward and reverse reads respectively. 
if PAIRED:
	SAMPLES = list(set([re.search("(^.+)_[12]$", x).group(1) for x in SAMPLES]))

if len(SAMPLES) == 0:
	sys.exit("ERROR: no samples in base directory!")

rule all:
	input:
		expand("bam/{sample}.bam.bai", sample = SAMPLES)

if PAIRED:
	rule trim:
		input:
			read1 = BASEDIR + "/{sample}_1" + EXTENSION,
			read2 = BASEDIR + "/{sample}_2" + EXTENSION
		output:
			out1 = "trimmed/{sample}_1.fastq",
			out2 = "trimmed/{sample}_2.fastq"
		params:
			adapter_file = ADAPTERS
		threads: 1
		log:
			"logs/{sample}_trim.log"
		shell:
			"bbduk.sh -Xmx4g -ignorebadquality in1={input.read1} in2={input.read2} out1={output.out1} out2={output.out2} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 > {log} 2>&1"

	rule bwa_mapping:
		input:
			trimmed1 = "trimmed/{sample}_1.fastq",
			trimmed2 = "trimmed/{sample}_2.fastq"
		output:
			"bwa_out/{sample}.sam"
		params:
			name = "{sample}",
			idx = BWA_INDEX
		log:
			"logs/{sample}_map.log"
		threads: 8
		shell: "bwa mem -t {threads} -R '@RG\tID:{params.name}\tSM:{params.name}' {params.idx} {input.trimmed1} {input.trimmed2} > {output}"

else:
	rule trim:
		input:
			BASEDIR + "/{sample}" + EXTENSION
		output:
			"trimmed/{sample}.fastq"
		params:
			adapter_file = ADAPTERS
		threads: 1
		log:
			"logs/{sample}_trim.log"
		shell:
			"bbduk.sh -Xmx4g -ignorebadquality in={input} out={output} ref={params.adapter_file} ktrim=r overwrite=true k=23 maq=20 mink=11 hdist=1 > {log} 2>&1"

	rule bwa_mapping:
		input:
			"trimmed/{sample}.fastq"
		output:
			"bwa_out/{sample}.sam"
		params:
			name = "{sample}",
			idx = BWA_INDEX
		log:
			"logs/{sample}_map.log"
		threads: 8
		shell: "bwa mem -t {threads} -R '@RG\tID:{params.name}\tSM:{params.name}' {params.idx} {input} > {output}"

rule convert_to_bam:
	input:
		"bwa_out/{sample}.sam"
	output:
		"bam/{sample}.bam"
	shell:
		"samtools view -bh {input} | samtools sort -o {output} -"

rule index_bam:
	input:
		"bam/{sample}.bam"
	output:
		"bam/{sample}.bam.bai"
	shell:
		"samtools index {input}"
