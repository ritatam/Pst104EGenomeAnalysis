#!/bin/bash

hap=hapA
wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/funannotate/external_annotate
tmhmm_results=$wdir/Puccinia_striiformis_Pst104E_${hap}.tmhmm.txt #make sure this is proteome
ingenannot_effector_results_custom=$wdir/Puccinia_striiformis_Pst104E_${hap}.ingenannot_effector_predictor-noTM.custom.txt
output=$wdir/Puccinia_striiformis_Pst104E_${hap}.final.custom.txt

# make custom tmhmm annotation.
tmhmm_custom=$wdir/Puccinia_striiformis_Pst104E_${hap}.tmhmm.custom.txt
 echo -n > $tmhmm_custom
 while read -r row
 do
     TM_gene=$(echo $row | awk '{print $1}')
     numHel=$(echo $row | awk '{print $5}' | sed 's/PredHel=//')
     topology=$(echo $row | awk '{print $6}' | sed 's/Topology=//')
     echo -e "${TM_gene}\tnote\tTransmembrane_tmhmm:${numHel} (${topology})" >> $tmhmm_custom
 done < <(grep -v "PredHel=0" $tmhmm_results)

# then combine with ingenannot effector noTM annotation to produce final custom annotation
# first check if there are genes common to both tmhmm and ingenannot effector annotations;
# if there is, will need to combine them into one row in the final custom annotation
common_genes=$(comm -12 <(cut -f1 ${ingenannot_effector_results_custom} | sort) <(cut -f1 ${tmhmm_custom} | sort))

if [  ! -z $common_genes ]; then
    echo "The following genes are common to ${tmhmm_custom} and ${ingenannot_effector_results_custom}:"
    echo "$common_genes"
    exit 1  # Exit immediately
fi
cat ${tmhmm_custom} ${ingenannot_effector_results_custom} > $output
