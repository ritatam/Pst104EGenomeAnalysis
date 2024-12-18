#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/common-tools/ 

craq=/media/ssd/rita/softwares/CRAQ/bin/craq

for hap in {A,B}
do
	# duplex
	echo chr{1..18}$hap | xargs samtools view -bh ../duplex.v3.9.chr.bam -@8 | samtools sort -@8 -o duplex.v3.9.hap${hap}.bam
	samtools index -@8 duplex.v3.9.hap${hap}.bam
	# illumina
	echo chr{1..18}$hap | xargs samtools view -bh ../illumina.v3.9.chr.bam -@8 | samtools sort -@8 -o illumina.v3.9.hap${hap}.bam
	samtools index -@8 illumina.v3.9.hap${hap}.bam
	$craq -t 10 -g ../v3.9.hap${hap}.fasta -D hap$hap -ngs illumina.v3.9.hap${hap}.bam -x map-ont -sms duplex.v3.9.hap${hap}.bam
done

$craq -t 10 -g ../v3.9.chr.haplotype-paired.fasta -D chr -x map-ont -ngs ../illumina.v3.9.chr.bam -sms ../duplex.v3.9.chr.bam
