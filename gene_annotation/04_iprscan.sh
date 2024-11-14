#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/funannotate

hap=hapA

export FUNANNOTATE_DB=/media/ssd/rita/project/104e/gene_annotations/funannotate/funannotate_db
export PATH=$PATH:/home/groups/schwessinger/condaEnvs/funannotate/opt/gmes_linux_64_4
export EGGNOG_DATA_DIR=/home/groups/schwessinger/condaEnvs/funannotate/opt/eggnog-mapper/data

wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/funannotate/external_annotate
ref=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/data/v3.9.${hap}.masked.fasta # make sure it's repeat-masked, sorted and cleaned
input_gff=/media/ssd/rita/project/104e/gene_annotations/ingenannot_utr/utr_refine/Puccinia_striiformis_Pst104E_${hap}.genes.utrs.trna_added.renamed.alias_rm.gff3


ln -s $input_gff $wdir/Puccinia_striiformis_Pst104E_${hap}.gff3

funannotate util gff2prot \
	--gff3 $wdir/Puccinia_striiformis_Pst104E_${hap}.gff3 \
	--fasta $ref \
	--no_stop > $wdir/Puccinia_striiformis_Pst104E_${hap}.proteins.fa

funannotate iprscan -i $wdir/Puccinia_striiformis_Pst104E_${hap}.proteins.fa \
	-m local \
	-o $wdir/Puccinia_striiformis_Pst104E_${hap}.iprscan.5.64-96.0.local.xml \
	-c 16 \
	--iprscan_path /media/ssd/rita/softwares/interproscan-5.64-96.0/interproscan.sh
