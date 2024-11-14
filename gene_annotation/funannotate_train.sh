#!/bin/bash


source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/funannotate

export FUNANNOTATE_DB=/media/ssd/rita/project/104e/gene_annotations/funannotate/funannotate_db
export GENEMARK_PATH=gmes/home/groups/schwessinger/condaEnvs/funannotate/opt/gmes_linux_64_4
export PATH=$PATH:/home/groups/schwessinger/condaEnvs/funannotate/opt/gmes_linux_64_4
export EGGNOG_DATA_DIR=/home/groups/schwessinger/condaEnvs/funannotate/opt/eggnog-mapper/data

wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/hapA
ref=/media/ssd/rita/project/104e/gene_annotations/funannotate/hapA/data/v3.9.hapA.masked.fasta # make sure it's repeat-masked, sorted and cleaned

cd $wdir
funannotate train -i $ref -o $wdir/funannotate \
	--left data/Pst104E_merged.sr.paired-mapped.hapA.r1.fastq.gz \
	--right data/Pst104E_merged.sr.paired-mapped.hapA.r2.fastq.gz \
        --trinity data/concatenated.espresso.stringtie2.trinity.hapA.fasta \
	--stranded no \
        --jaccard_clip \
        --species "Puccinia striiformis" \
        --isolate Pst104E_hapA \
        --no_trimmomatic \
        --cpus 14 \
        --memory 100G \
        --pasa_db mysql 
