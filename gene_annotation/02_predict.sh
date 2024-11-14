#!/bin/bash


source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/funannotate

export FUNANNOTATE_DB=/media/ssd/rita/project/104e/gene_annotations/funannotate/funannotate_db
export PATH=$PATH:/home/groups/schwessinger/condaEnvs/funannotate/opt/gmes_linux_64_4
export EGGNOG_DATA_DIR=/home/groups/schwessinger/condaEnvs/funannotate/opt/eggnog-mapper/data

wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/hapA
ref=/media/ssd/rita/project/104e/gene_annotations/funannotate/hapA/data/v3.9.hapA.masked.fasta # make sure it's repeat-masked, sorted and cleaned


funannotate predict -i $ref -o $wdir/funannotate \
	--species "Puccinia striiformis" \
	--isolate Pst104E_hapA \
	--name Pst104E137 \
	--transcript_evidence $wdir/data/concatenated.espresso.stringtie2.trinity.hapA.fasta \
	--rna_bam $wdir/funannotate/training/funannotate_train.coordSorted.bam \
	--pasa_gff $wdir/funannotate/training/funannotate_train.pasa.gff3 \
	--stringtie $wdir/funannotate/training/funannotate_train.stringtie.gtf \
	--other_gff $wdir/data/codingquarry-pm.hapA.gff3:10 \
        --transcript_alignments ${wdir}/data/stringtie_merged.espresso-stringtie2-stringtie.transcripts.hapA.gff3 \
	--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta $wdir/data/Pst_104E_v13_ph_ctg.protein.fa \
	--weights augustus:4 hiq:6 genemark:1 pasa:10 codingquarry:0 snap:1 glimmerhmm:1 proteins:6 transcripts:6 \
	--optimize_augustus \
	--repeats2evm \
	--ploidy 1 \
	--cpus 30

