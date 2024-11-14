# signalp6
signalp6 --fasta Puccinia_striiformis_Pst104E_hapA.proteins.fa --organism eukarya --output_dir signlap6/ --format txt --mode slow

# tmhmm
# input proteome for transmembrane domain annotation
tmhmm external_annotate/Puccinia_striiformis_Pst104E_hapB.proteins.fa -short > proteome.tmhmm.hapB.txt
# input signalp6 SP-cleaved outputs for filtering false positives
tmhmm signalP6.hapA.faa  -short > signalP6.tmhmm.hapA.txt

# phobius - ran on web server, once on proteome, once only on signalP6 SP-cleaved outputs

