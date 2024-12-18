#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/merqury

input_fastq=/media/nvme/rita/data/nanopore/104e/dorado/Pst104E_dorado_duplex.q30.l10k.htcrop75.yacrd_split.fastq
outdir=/media/ssd/rita/project/104e/merqury/v3.9_hapA
ref=/media/ssd/rita/project/104e/merqury/v3.9_chr/v3.9.chr.haplotype-paired.fasta
mkdir -p $outdir
cd $outdir

# 1. reads fastq to fasta
reads_fasta=${input_fastq%.fastq}.fasta
#seqtk seq -A ${input_fastq} > ${reads_fasta}


#conda deactivate
# 2. meryl count kmers
export PATH=/media/nvme/rita/softwares/meryl-1.4/build/bin:$PATH
meryl count k=31 ${reads_fasta} output duplex.meryl threads=10 memory=32g

# 3. merqury
export MERQURY=/home/groups/schwessinger/condaEnvs/merqury/share/merqury
$MERQURY/merqury.sh duplex.meryl $ref  duplex_merqury
