#!/bin/bash

# run verkko assembly + phase assembly with verkkohic pipeline.

source /opt/conda/etc/profile.d/conda.sh

### SETTINGS ###
workdir=/media/ssd/rita/project/104e/verkko_asm
duplex=${workdir}/hifi/Pst104E_duplex.q30.l10k.htcrop75.yacrd_split.fastq.gz
simplex=${workdir}/ont/Pst104E_simplex.q10.l100k.htcrop100.fastq.gz
unphased_out=${workdir}/unphased_asm
################

echo Duplex is $duplex
echo Simplex is $simplex

# run verkko
conda activate /home/groups/schwessinger/condaEnvs/verkko ;
mkdir -p ${unphased_out} ;
verkko -d ${unphased_out} --hifi ${duplex} --nano ${simplex} --threads 16 --local-memory 56 --mbg /media/ssd/rita/softwares/MBG/bin/MBG ;

# run verkkohic

conda deactivate
conda activate /home/groups/schwessinger/condaEnvs/verkko_tmp

mkdir -p ${gfase_out}
cd ${workdir}
gfase_wrapper=/media/ssd/rita/softwares/verkkohic/gfase_wrapper.sh
export VERKKO=/media/ssd/rita/softwares/verkko
export GFASE=/media/ssd/rita/softwares/GFAse/build
export PATH="/media/ssd/rita/softwares/MBG/bin:$PATH"
export PATH="/media/ssd/rita/softwares/verkko/bin:$PATH"

bash $gfase_wrapper unphased_asm gfase_asm `pwd`