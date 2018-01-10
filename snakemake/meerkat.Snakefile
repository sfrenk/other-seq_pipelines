import glob
import re
import sys

############################    DESCRIPTION    ##############################

# Run Meerkat on bam files
# NOTE: use setup_dir.sh -d <directory_name> -p gatk to set up the working directory

############################    PARAMETERS    ##############################

# Input parameters
BASEDIR = "/nas/longleaf/home/sfrenk/pine/peter_carlton/paired_end/bam/"

###############################################################################

SAMPLE_FILES = glob.glob(BASEDIR + "*" + ".bam")
SAMPLES = [ re.search(BASEDIR + "/?([^/]+)" + ".bam", x).group(1) for x in SAMPLE_FILES ]

if len(SAMPLES) == 0:
	sys.exit("ERROR: no samples in base directory!")

rule all:
	input:
		expand("{sample}.variants", sample = SAMPLES)

rule pre_process:
	input:
		BASEDIR + "{sample}.bam"
	output:
		BASEDIR + "{sample}.isinfo"
	shell:
		"perl /nas/longleaf/home/sfrenk/local/Meerkat/scripts/pre_process.pl -s 20 -k 1500 -q 15 -l 0 -b {input}"

rule run_meerkat:
	input:
		bamfile = BASEDIR + "{sample}.bam",
		infofile = BASEDIR + "{sample}.isinfo"
	output:
		"{sample}.variants"
	threads: 8
	run:
		shell("perl /nas/longleaf/home/sfrenk/local/Meerkat/scripts/meerkat.pl -s 20 -d 5 -p 3 -o 1 -m 0 -l 0 -t 8 -F /proj/ahmedlab/steve/Software/Meerkat/WS251/fasta -b {input.bamfile}")
		shell("perl /nas/longleaf/home/sfrenk/local/Meerkat/scripts/mechanism.pl -R /proj/ahmedlab/steve/Software/Meerkat/WS251/ce11_rmsk.txt -b {input.bamfile}")
