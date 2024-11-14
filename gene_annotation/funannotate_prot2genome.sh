#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/funannotate

hap=hapA

export FUNANNOTATE_DB=/media/ssd/rita/project/104e/gene_annotations/funannotate/funannotate_db
wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/funannotate/external_annotate
ref=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/data/v3.9.${hap}.masked.fasta # make sure it's repeat-masked, sorted and cleaned

funannotate util prot2genome \
	-g $ref \
	-p $wdir/Puccinia_striiformis_Pst104E_${hap}.proteins.fa \
	-o $wdir/Puccinia_striiformis_Pst104E_${hap}.p2g.gff3 \
	--cpus 24 \
	--logfile ${hap}-p2g.log

