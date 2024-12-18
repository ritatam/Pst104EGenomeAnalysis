# raw hic reads must be in a input folder, inside another subdirectory named by the isolate or else like so work_dir/HiC_input/Pst104E/HiC_R*.fastq

# generate list of restriction fragments from digested genome.
HICPRO_UTIL=/media/nvme/rita/softwares/HiC-Pro/bin/utils
$HICPRO_UTILS/digest_genome.py -r ^GATC,G^ANTC,T^TAA,C^TNAG -o assembly.v3.9.digested.bed assembly.v3.9.fasta

# generate bowtie2 index file for the genome.
bowtie2-build assembly.v3.9.fasta bowtie2_index/asssembly.v3.9

# generate a table file of chr sizes.
samtools faidx assembly.v3.9.fasta
cut -f1,2 assembly.v3.9.fasta.fai > assembly.v3.9.chrsize

# edits on config-hicpro.txt:
# MIN_MAPQ = 20
# BOWTIE2_IDX_PATH = 104e_verkko_v2/hicpro/bowtie2_index
# GENOME_FRAGMENT = 104e_verkko_v2/hicpro/assembly.v3.9.digested.bed
# LIGATION_SITE = GATCGATC,GANTANTC,TTATAA,CTNATNAG
# REFERENCE_GENOME = assembly.v3.9   #This must be the prefix used in BOWTIE2_IDX_PATH
# GENOME_SIZE = 104e_verkko_v2/hicpro/assembly.v3.9.chrsize

HiC-Pro -i HiC_input -o output -c config-hicpro.txt
