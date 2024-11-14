#!/bin/bash

# This script performs funannotate gene gff3 processing to meet formatting requirements by ingenannot.

source /opt/conda/etc/profile.d/conda.sh

hap=hapA
funannotategff=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/funannotate/predict_results/Puccinia_striiformis_Pst104E_${hap}.gff3
genegff=../../ingenannot_rescue_effectors/Puccinia_striiformis_Pst104E_${hap}.trna_commented.rescued_effectors.gff3
isoform_rank=../isoform_rank/isoforms_${hap}.alternatives.gff
output=Puccinia_striiformis_Pst104E_${hap}.genes.utrs.gff

# search for tRNA feature and writes out the corresponding gene id
echo "writing out tRNA from $funannotategff"
grep tRNA $funannotategff | cut -f9 | sed -n "s/ID=\([^;]*\)-.*/\1/p" > tRNA_genes.${hap}.list

# search for tRNA genes in ingenannot's isoform_rank.alternatives.gff output
while read -r gene_id
do
        if grep -q $gene_id $isoform_rank
        then
                echo -e "tRNA found for gene ID $gene_id in isoform_rank.alternatives.gff, please edit it out\n"
                exit 1
        fi
done < tRNA_genes.${hap}.list
echo -e "No isoform entries for tRNA genes found in ${isoform_rank}, continuing pipeline \n"


# use the tmp file to add UTRs using isoform rank as preferred isoform
conda activate /home/groups/schwessinger/condaEnvs/ingenannot
echo -e "adding UTRs using isoform rank as preferred isoform"
ingenannot -v 2 utr_refine $genegff $isoform_rank $output --erase --utr_mode rank

echo -e "concatenating $output with tRNA_entries.${hap}.gff3"
grep "^#chr" $genegff | sed "s/#chr/chr/g" > tRNA_entries.${hap}.gff3
cat $output tRNA_entries.${hap}.gff3 > ${output%.gff}.trna_added.gff
