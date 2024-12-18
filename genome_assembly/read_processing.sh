### duplex reads processing ###
raw_duplex=Pst104E_duplex.fastq
filt_duplex=Pst104E_duplex.q30.l10k.htcrop75.fastq.gz
paf=Pst104E_duplex.q30.l10k.htcrop75.ava.paf
chimera_split_duplex=Pst104E_duplex.q30.l10k.htcrop75.yacrd_split.fastq

# filter and hard-trim duplex 
NanoFilt -q 30 -l 10000 --headcrop 75 --tailcrop 75 $raw_duplex | gzip > $filt_duplex
# all-vs-all duplex alignment, then fpa to drop mapping with length shorter than 1000 bp.
minimap2 -x ava-ont -t 48 $filt_duplex $filt_duplex | fpa drop -l 1000 > $paf
# conservatively split chimera where ava coverage is 0 based on yacrd author's recommendation (https://github.com/natir/yacrd/issues/45)
yacrd -i $paf -o report.yacrd split -i $filt_duplex -o $chimera_split_duplex


### simplex reads processing ###
raw_simplex=Pst104E_simplex.fastq
filt_simplex=Pst104E_simplex.q10.l40k.htcrop100.fastq.gz
NanoFilt -q 10 -l 40000 --headcrop 100 --tailcrop 100 $raw_simplex | gzip > $filt_simplex


