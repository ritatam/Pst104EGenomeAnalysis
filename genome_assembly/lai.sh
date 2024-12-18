#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/ltr_retriever
gt=/media/ssd/rita/softwares/LTR_retriever/gt-1.6.2-Linux_x86_64-64bit-complete/bin/gt
ltr_finder_parallel=/media/ssd/rita/softwares/LTR_retriever/LTR_FINDER_parallel/LTR_FINDER_parallel
ltr_retriever=/media/ssd/rita/softwares/LTR_retriever/LTR_retriever
lai=/media/ssd/rita/softwares/LTR_retriever/LAI

for hap in {hapA,hapB}
do
	genome=/media/ssd/rita/project/104e/LAI/v3.9.${hap}.fasta

	$gt suffixerator -db $genome -indexname $genome -tis -suf -lcp -des -ssp -sds -dna

	$gt ltrharvest -index $genome -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > ${genome}.harvest.scn

	$ltr_finder_parallel -seq $genome -threads 10 -harvest_out -size 1000000 -time 300
	
	cat $genome.harvest.scn $genome.finder.combine.scn > $genome.rawLTR.scn
        
	$ltr_retriever -genome $genome -inharvest $genome.rawLTR.scn -threads 24

done

