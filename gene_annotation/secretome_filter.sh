#!/bin/bash

# filter signalP 6.0 secretomes based on TM domains detected by BOTH phobius and tmhmm to ensure it's biologically accurate.
# it also filters the ingenannot effector by applying the same filter, if any has TM detected by both. this should be unlikely.

source /opt/conda/etc/profile.d/conda.sh

hap=hapA
wdir=/media/ssd/rita/project/104e/gene_annotations/funannotate/${hap}/funannotate/external_annotate/secretome_pred

signalp6_results=$wdir/../Puccinia_striiformis_Pst104E_${hap}.signalp6.txt
tmhmm_results=$wdir/tmhmm_results.signalp6aa.${hap}.txt #signalp6 entries only
phobius_results=$wdir/phobius_results.signalp6aa.${hap}.txt #signalp6 entries only
ingenannot_effectors=$wdir/../Puccinia_striiformis_Pst104E_${hap}.ingenannot_effector_predictor.custom.txt

# extract gene IDs of SPs that have TM detected by both TMHMM and phobius
mkdir -p SP-TM
tmhmm_TM=$wdir/SP-TM/tmhmm-TM-SPs.${hap}.names
phobius_TM=$wdir/SP-TM/phobius-TM-SPs.${hap}.names
TM=$wdir/SP-TM/TM-SPs.${hap}.names
awk '$5 != "PredHel=0" {print $1}' $tmhmm_results > $tmhmm_TM
awk '$2 != "0" {print $1}' <(tail -n+2 $phobius_results) > $phobius_TM

# intersect both TM lists. this will list genes with TM domain detected by both phobius and TMHMM.
comm -12 <(sort $tmhmm_TM) <(sort $phobius_TM)  > $TM

# extract signalP SP gene names and remove those that appear in $TM file.
# note "noTM" here only means no TM confirmed, not "truly biologically absent" - one of the tools might have detected TM but the other didn't, which doesn't satisfy the criteria we set up here
noTM=$wdir/noTM-SPs.${hap}.names
comm -3 <(cut -f1 $tmhmm_results | sort) <(sort $TM) > $noTM

# use the TM list to filter SPs from signalP6 and write to the signalP6-TM txt for funannotate
signalp6_noTM_results=$wdir/../Puccinia_striiformis_Pst104E_${hap}.signalp6-noTM.txt
head -n2 $signalp6_results > $signalp6_noTM_results
while read -r gene_id; do grep $gene_id <(tail -n+3 $signalp6_results) >> $signalp6_noTM_results; done<$noTM

# also remove effectors that have TM detected by both tools
ingenannot_noTM_results=$wdir/../Puccinia_striiformis_Pst104E_${hap}.ingenannot_effector_predictor-noTM.custom.txt
echo -n > $ingenannot_noTM_results
awk 'FNR==NR { names[$1]; next } !($1 in names)' $TM $ingenannot_effectors > $ingenannot_noTM_results
