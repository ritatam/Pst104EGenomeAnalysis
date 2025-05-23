#!/bin/bash
set -euox
source /opt/conda/etc/profile.d/conda.sh
conda activate /home/groups/schwessinger/condaEnvs/3d-dna
export JAVA_TOOL_OPTIONS="-Xmx150g"
export wdir=/media/nvme/rita/project/104e_verkko_v2/juicer

cd $wdir/3d-dna
bash /home/groups/schwessinger/condaEnvs/3d-dna/3d-dna/run-asm-pipeline.sh --mode diploid --rounds -1 --build-gapped-map --sort-output $wdir/references/assembly.v3.9.fasta $wdir/aligned/merged_nodups.txt

# take the heatmap generated by the first iteration (.0) for review in juicebox. 

# post review
bash /home/groups/schwessinger/condaEnvs/3d-dna/3d-dna/run-asm-pipeline-post-review.sh --review assembly.v3.9.0.review_v2.assembly $wdir/references/assembly.v3.9.fasta $wdir/aligned/merged_nodups.txt