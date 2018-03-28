import glob
import re
import sys

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = ""
EXTENSION = ""
PAIRED = False

# Trimming parameters
ADAPTERS = "~/proj/seq/bbmap/adapters.fa"

# Mapping parameters
BOWTIE_INDEX = "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/bowtie/genome"
MULTI_FLAG = "-M 1"
LEFT_TRIM = 0
RIGHT_TRIM = 0

# Coverage parameters
# Length by which to extend reads (leave blank for paired end)
EXT_LENGTH = ""
# Name for merged bg file
BG_NAME = "bg/merged.bg"

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
		BG_NAME

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

	rule bowtie_mapping:
		input:
			trimmed1 = "trimmed/{sample}_1.fastq",
			trimmed2 = "trimmed/{sample}_2.fastq"
		output:
			"bowtie_out/{sample}.sam"
		params:
			idx = BOWTIE_INDEX,
			multi_flag = MULTI_FLAG,
			left_trim = LEFT_TRIM,
			right_trim = RIGHT_TRIM
		log:
			"logs/{sample}_map.log"
		threads: 8
		shell: "bowtie -p {threads} -S {params.multi_flag} --chunkmbs 200 --best --strata -5 {params.left_trim} -3 {params.right_trim} {params.idx} --minins 100 --maxins 1000 -1 {input.trimmed1} -2 {input.trimmed2} {output} > {log} 2>&1"


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

	rule bowtie_mapping:
		input:
			"trimmed/{sample}.fastq"
		output:
			"bowtie_out/{sample}.sam"
		params:
			idx = BOWTIE_INDEX,
			multi_flag = MULTI_FLAG,
			left_trim = LEFT_TRIM,
			right_trim = RIGHT_TRIM
		log:
			"logs/{sample}_map.log"
		threads: 8
		shell: "bowtie -p {threads} -S {params.multi_flag} --chunkmbs 200 --best --strata -5 {params.left_trim} -3 {params.right_trim} {params.idx} {input} {output} > {log} 2>&1"

rule convert_to_bam:
	input:
		"bowtie_out/{sample}.sam"
	output:
		"bam/{sample}.bam"
	shell:
		"samtools view -bh -F 4 {input} | samtools sort -o {output} -"

rule index_bam:
	input:
		"bam/{sample}.bam"
	output:
		"bam/{sample}.bam.bai"
	shell:
		"samtools index {input}"

rule make_bg:
	input:
		bamfile = "bam/{sample}.bam",
		bamidx = "bam/{sample}.bam.bai"
	output:
		"bg/{sample}.bg"
	params:
		ext_length = EXT_LENGTH
	threads:
		8
	log:
		"logs/{sample}_bg.log"
	shell:
		"module purge; "
		"module load deeptools/2.5.4; "
		"bamCoverage -of bedgraph -bs 200 -b {input.bamfile} -o {output} -p {threads} --normalizeTo1x 100258171 --extendReads {params.ext_length} >> {log} 2>&1"

rule merge_bg:
	input:
		expand("bg/{sample}.bg", sample = SAMPLES)
	output:
		BG_NAME
	params:
		name = expand("{sample}", sample = SAMPLES)
	log:
		"logs/merge_bg.log"
	shell:
		"module load bedtools; "
		"bedtools unionbedg -i {input} -header -names {params.name} > {output} 2> {log}"
