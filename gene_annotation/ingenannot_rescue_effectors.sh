#!/bin/bash

hap=hapA

source /opt/conda/etc/profile.d/conda.sh

# validate data
conda activate /home/groups/schwessinger/condaEnvs/ingenannot
ingenannot -v 2 validate Puccinia_striiformis_Pst104E_${hap}.trna_commented.gff3
ingenannot -v 2 validate stringtie_merged.espresso-stringtie2-stringtie.transcripts.${hap}.gtf

# sort and index transcript
conda deactivate
conda activate /home/groups/schwessinger/condaEnvs/transcript-anno
sort -k1,1 -k4g,4 stringtie_merged.espresso-stringtie2-stringtie.transcripts.${hap}.gtf > stringtie_merged.espresso-stringtie2-stringtie.transcripts.${hap}.sorted.gtf
bgzip stringtie_merged.espresso-stringtie2-stringtie.transcripts.${hap}.sorted.gtf
tabix stringtie_merged.espresso-stringtie2-stringtie.transcripts.${hap}.sorted.gtf.gz

# run rescue effector on unused transcripts 
conda deactivate
conda activate /home/groups/schwessinger/condaEnvs/ingenannot
export PATH=$PATH:/media/ssd/rita/softwares/signalp-4.1
export PATH=$PATH:/media/ssd/rita/softwares/targetp-1.1 #remember to configure targetp file properly
export PATH=$PATH:/media/ssd/rita/softwares/tmhmm-2.0c/bin # all tmhmm perl scripts must have shebang edited to the correct perl path

ingenannot -v 2 rescue_effectors Puccinia_striiformis_Pst104E_${hap}.trna_commented.gff3 stringtie_merged.espresso-stringtie2-stringtie.transcripts.${hap}.sorted.gtf.gz ../../../assembly_versions/v3.9_gapfill/v3.9.${hap}.fasta --effectorp /media/ssd/rita/softwares/EffectorP-2.0/Scripts/EffectorP.py 

# concatenate effectors.gff3 and funannotate-predicted gene gff3.
# use funannotate util for proper concatenation and numbering.
# note this will remove the masked tRNA entries so will need to add them back afterwards.
cat Puccinia_striiformis_Pst104E_hapA.trna_commented.gff3 effectors.gff3 > ../Puccinia_striiformis_Pst104E_hapA.trna_commented.rescued_effectors.gff3
cat Puccinia_striiformis_Pst104E_hapB.trna_commented.gff3 effectors.gff3 > ../Puccinia_striiformis_Pst104E_hapB.trna_commented.rescued_effectors.gff3

funannotate util gff-rename --gff3 Puccinia_striiformis_Pst104E_hapA.trna_commented.rescued_effectors.gff3 --fasta ../../assembly_versions/v3.9_gapfill/v3.9.hapA.fasta --locus_tag Pst104E137 --out Puccinia_striiformis_Pst104E_hapA.trna_commented.rescued_effectors.renamed.gff3

funannotate util gff-rename --gff3 Puccinia_striiformis_Pst104E_hapB.trna_commented.rescued_effectors.gff3 --fasta ../../assembly_versions/v3.9_gapfill/v3.9.hapB.fasta --locus_tag Pst104E137 --out Puccinia_striiformis_Pst104E_hapB.trna_commented.rescued_effectors.renamed.gff3
