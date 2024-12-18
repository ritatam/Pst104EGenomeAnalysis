**clean up mtDNA contigs post blast**

First concatenate all blast outputs into one file

```
cat *.blast > ../all.blast
```

First, get names of chunks/contigs that have keyword "mitochon" in their blast output
```
grep mitochon all.blast > mtDNA.blast
cut -f1 mtDNA.blast | sort | uniq > mtDNA_contigs.txt
wc -l mtDNA_contigs.txt
```

**I strong recommend users to examine qcov and alignment length against your target species/genus mtDNA records in nt database (eg Pucciniales), to double check if there's false positive. For example a mtDNA hit with very short qcov could indicate shared sequences between gDNA and mtDNA e.g. cytochrome b gene.**

Once inspected and confirmed, save mtDNA contigs in mtDNA.fasta to keep record. 
```
seqtk subseq assembly.fasta mtDNA_contigs.txt > mtDNA.fasta
```

use bbmap to remove mtDNA contigs from assembly
```
filterbyname.sh in=gfase_asm/assembly.fasta out=assembly.mtDNArm.fasta names=mtDNA_contigs.txt
```

**clean up contaminant contigs**

check all species names 
`cut -f12 all.blast | sort | uniq | less`

Printed:
```
Austropuccinia psidii
Puccinia aucta
Puccinia cf. psidii AE-2014
Puccinia coronata f. sp. avenae
Puccinia graminis f. sp. tritici
Puccinia graminis f. sp. tritici CRL 75-36-700-3
Puccinia hordei
Puccinia poarum
Puccinia striiformis f. sp. hordei
Puccinia striiformis f. sp. tritici
Puccinia triticina
uncultured fungus
Uredo kriegeriana
Uromyces dactylidis
```
All are fungal plant pathogen (uncultured fungus).
grep species keyword from all.blast. If you see rDNA hits from other closely related species it's pretty normal

**In case there's contaminant assembled into contigs:**
BLAST tells which part of the contig is blasting with a contaminant (or non rust DNA) and gives bp positions. In this case one may want to map the reads back to the genome assembly, and visualise the alignments at the suspected contaminant location to double see if there's alignment cutoff & check boundaries. 