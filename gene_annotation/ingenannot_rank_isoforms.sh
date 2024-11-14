#!/bin/bash

hap=hapA
transcriptdir=/media/ssd/rita/project/104e/gene_annotations/transcripts
bamdir=/media/longterm/rita/104e_nanopore_cDNA/mapped/Pst104Ev3.9/mapped/merged/${hap}
source /opt/conda/etc/profile.d/conda.sh

# merge the long-read alignment bam files for coverage and junction support
conda activate /home/groups/schwessinger/condaEnvs/common-tools
merged_bam=Pst104E_all_merged.ont.mapped.${hap}.bam
sorted_bam=Pst104E_all_merged.ont.mapped.${hap}.sorted.bam
echo -e "samtools merging long-read alignment bam files"
samtools merge -o $merged_bam $(ls ${bamdir}/*.bam) -@16
samtools sort -@16 -O BAM $merged_bam -o $sorted_bam
samtools index -@16 $sorted_bam
conda deactivate

mkdir -p isoform_rank
# rank RNA-seq and ONT cDNA transcript annotations
echo -e "rank RNA-seq and ONT cDNA transcript annotations"
conda activate /home/groups/schwessinger/condaEnvs/ingenannot
ingenannot -v 2 isoform_ranking stringtie_merged.espresso-stringtie2-stringtie.transcripts.${hap}.gtf -b $sorted_bam  --alt_threshold 0.1 --rescue --prefix isoform_rank/isoforms_${hap}
