conda activate /home/groups/schwessinger/condaEnvs/ingenannot
export PATH=$PATH:/media/ssd/rita/softwares/signalp-4.1
export PATH=$PATH:/media/ssd/rita/softwares/targetp-1.1
export PATH=$PATH:/media/ssd/rita/softwares/tmhmm-2.0c/bin 
ingenannot -v 2 effector_predictor Puccinia_striiformis_Pst104E_hapA.proteins.fa --effectorp /media/ssd/rita/softwares/EffectorP-2.0/Scripts/EffectorP.py
ingenannot -v 2 effector_predictor Puccinia_striiformis_Pst104E_hapB.proteins.fa --effectorp /media/ssd/rita/softwares/EffectorP-2.0/Scripts/EffectorP.py
# format it as custom annotation txt file for funannotate
tail -n +4 effectors.txt | while IFS="," read -ra genes; do for gene in ${genes[0]} ; do echo -e "$gene\tnote\tingenannot:predicted_effector" ; done; done > Puccinia_striiformis_Pst104E_hapA.ingenannot_effector_predictor.custom.txt
tail -n +4 effectors.txt | while IFS="," read -ra genes; do for gene in ${genes[0]} ; do echo -e "$gene\tnote\tingenannot:predicted_effector" ; done; done > Puccinia_striiformis_Pst104E_hapB.ingenannot_effector_predictor.custom.txt
