#!/bin/bash

hap=hapA

source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/funannotate

export FUNANNOTATE_DB=/media/ssd/rita/project/104e/gene_annotations/funannotate/funannotate_db
export GENEMARK_PATH=gmes/home/groups/schwessinger/condaEnvs/funannotate/opt/gmes_linux_64_4
export PATH=$PATH:/home/groups/schwessinger/condaEnvs/funannotate/opt/gmes_linux_64_4
export EGGNOG_DATA_DIR=/home/groups/schwessinger/condaEnvs/funannotate/opt/eggnog-mapper/data

wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/hapA/funannotate
ref=/media/ssd/rita/project/104e/gene_annotations/funannotate/hapA/data/v3.9.hapA.masked.fasta # make sure it's repeat-masked, sorted and cleaned

funannotate annotate \
	--gff $wdir/external_annotate/Puccinia_striiformis_Pst104E_${hap}.gff3 \
	--fasta $ref \
	--species "Puccinia striiformis" \
	--out $wdir \
	--antismash $wdir/external_annotate/Puccinia_striiformis_Pst104E_${hap}.antismash.gbk \
	--iprscan $wdir/external_annotate/Puccinia_striiformis_Pst104E_${hap}.iprscan.5.64-96.0.local.xml \
	--phobius $wdir/external_annotate/Puccinia_striiformis_Pst104E_${hap}.phobius.txt \
	--signalp $wdir/external_annotate/Puccinia_striiformis_Pst104E_${hap}.signalp6-noTM.txt \
	--p2g $wdir/external_annotate/Puccinia_striiformis_Pst104E_${hap}.p2g.gff3 \
	--annotations $wdir/external_annotate/Puccinia_striiformis_Pst104E_${hap}.final.custom.txt \
	--isolate Pst104E_hapA \
	--sbt $wdir/../../Puccinia_striiformis_Pst104E.sbt \
	--header_length 20 \
	--busco_db basidiomycota \
	--cpus 30 
