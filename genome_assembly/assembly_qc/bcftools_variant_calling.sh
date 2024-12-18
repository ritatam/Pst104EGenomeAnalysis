#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh

conda activate /home/groups/schwessinger/condaEnvs/common-tools
# filter for only primary alignments
samtools view -@24 -b -F 0x800 -F 0x100 ../duplex.v3.9.chr.bam -o duplex.v3.9.chr.primary-aln.bam
samtools index -@24 duplex.v3.9.chr.primary-aln.bam
conda deactivate
conda activate /home/groups/schwessinger/condaEnvs/SNP
bcftools mpileup --threads 24 -Ou -f ../chromosomes/v3.9.chr.haplotype-paired.fasta duplex.v3.9.chr.primary-aln.bam | bcftools call -mv --skip-variants indels -Ob -o duplex.v3.9.chr.SNPs.bcf

bcftools view -i 'QUAL>=20' duplex.v3.9.chr.SNPs.bcf > duplex.v3.9.chr.SNPs.q20.vcf


