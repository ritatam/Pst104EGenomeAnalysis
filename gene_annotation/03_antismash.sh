#!/bin/bash

hap=hapA

source /opt/conda/etc/profile.d/conda.sh
export FUNANNOTATE_DB=/media/ssd/rita/project/104e/gene_annotations/funannotate/funannotate_db

wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/funannotate/external_annotate
ref=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/data/v3.9.${hap}.masked.fasta # make sure it's repeat-masked, sorted and cleaned


wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/funannotate/external_annotate
ref=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/data/v3.9.${hap}.masked.fasta # make sure it's repeat-masked, sorted and cleaned
input_gff=/media/ssd/rita/project/104e/gene_annotations/ingenannot_utr/utr_refine/Puccinia_striiformis_Pst104E_hapA.genes.utrs.trna_added.renamed.alias_rm.gff3


ln -s $input_gff $wdir/Puccinia_striiformis_Pst104E_${hap}.gff3

conda activate /home/groups/schwessinger/condaEnvs/funannotate
mkdir -p $wdir/antismash
# convert gff to tbl
funannotate util gff2tbl \
	-g $wdir/Puccinia_striiformis_Pst104E_${hap}.gff3 \
	-f $ref > $wdir/antismash/Puccinia_striiformis_Pst104E_${hap}.tbl

# convert tbl to gbk
funannotate util tbl2gbk \
	--tbl $wdir/antismash/Puccinia_striiformis_Pst104E_${hap}.tbl \
	-f $ref \
	--species "Puccinia striiformis" \
	--isolate Pst104E \
	-o $wdir/antismash/Puccinia_striiformis_Pst104E_${hap}

conda deactivate

# antismash
conda activate /home/groups/schwessinger/condaEnvs/antismash
antismash --taxon fungi --output-dir $wdir/antismash/antismash --output-basename Puccinia_striiformis_Pst104E_${hap} --smcog-trees --cb-knownclusters --asf --cb-subclusters $wdir/antismash/Puccinia_striiformis_Pst104E_${hap}.gbk

echo -e "\ndone!"
cp $wdir/antismash/antismash/Puccinia_striiformis_Pst104E_${hap}.gbk $wdir/Puccinia_striiformis_Pst104E_${hap}.antismash.gbk
