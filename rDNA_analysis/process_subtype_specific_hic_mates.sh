#!/bin/bash

kmer_fasta=$1    # rDNA subtype-specific or unspecific k-mers
subtype=$2       # rDNA subtype name (will be used for output dir)
genome_fasta=$3  # reference for hic mates to map against
r1=$4            # hic r1
r2=$5            # hic r2
cpus=$6          # number of cpus to use for parallel run

source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/common-tools

# index genome
#bwa-mem2 index -p $genome_fasta $genome_fasta
# index kmer
awk '{ if ($0 ~ /^>/) { if (NR > 1) printf("\n"); printf("%s\t", substr($0, 2)); } else printf("%s", $0) } END { printf("\n") }' $kmer_fasta > $kmer_fasta.idx

process_kmer() {
    kmer_id=$1
    kmer=$2
    subtype=$3
    genome_fasta=$4
    r1=$5
    r2=$6

    mkdir -p $subtype
    outdir=$subtype/kmer_${kmer_id}
    mkdir -p $outdir
    r1_out=kmer_${kmer_id}.R1
    r2_out=kmer_${kmer_id}.R2

    # extract r1/r2 reads containing the specified unique kmer
    grep $kmer $r1 > ${outdir}/${r1_out}.fastq -B1 -A2
    grep $kmer $r2 > ${outdir}/${r2_out}.fastq -B1 -A2
    # grep and format identifiers of extracted reads associated with the kmer
    grep '^@' ${outdir}/${r1_out}.fastq | cut -d ' ' -f1 | sed 's/^@//' > ${outdir}/kmer_${kmer_id}.R1.readID.txt
    grep '^@' ${outdir}/${r2_out}.fastq | cut -d ' ' -f1 | sed 's/^@//' > ${outdir}/kmer_${kmer_id}.R2.readID.txt
    # extract mates
    seqtk subseq $r2 ${outdir}/kmer_${kmer_id}.R1.readID.txt > ${outdir}/${r2_out}.mates.fastq
    seqtk subseq $r1 ${outdir}/kmer_${kmer_id}.R2.readID.txt > ${outdir}/${r1_out}.mates.fastq
    # map mates to phased ref
    bwa-mem2 mem -t 16 $genome_fasta ${outdir}/${r2_out}.mates.fastq | samtools sort -@ 4 -o ${outdir}/kmer_${kmer_id}.R2.mates.bam
    bwa-mem2 mem -t 16 $genome_fasta ${outdir}/${r1_out}.mates.fastq | samtools sort -@ 4 -o ${outdir}/kmer_${kmer_id}.R1.mates.bam
    # convert bam to bed for mapping locations and MAPQ score
    # remove duplicated records as some reads might be tagged by consecutive k-mers multiple times
    bedtools bamtobed -i ${outdir}/kmer_${kmer_id}.R2.mates.bam | sort -k4,4 | uniq > ${outdir}/kmer_${kmer_id}.R2.mates.bed
    bedtools bamtobed -i ${outdir}/kmer_${kmer_id}.R1.mates.bam | sort -k4,4 | uniq > ${outdir}/kmer_${kmer_id}.R1.mates.bed
}

export -f process_kmer

# generate cmds and parallelise
awk -v subtype=$subtype -v genome=$genome_fasta -v r1=$r1 -v r2=$r2 \
    '{ print $1, $2, subtype, genome, r1, r2 }' $kmer_fasta.idx | parallel -j $cpus --colsep ' ' process_kmer