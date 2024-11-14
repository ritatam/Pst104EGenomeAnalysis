# order and renumber the protein-coding + tRNA + rescued effector genes predictions from funannotate and ingenannot, before proceeding to functional annotation.

#hapA
funannotate util gff-rename --gff3 Puccinia_striiformis_Pst104E_hapA.genes.utrs.trna_added.gff --fasta ../../../assembly_versions/v3.9_gapfill/v3.9.hapA.fasta --locus_tag Pst104E137 --out Puccinia_striiformis_Pst104E_hapA.genes.utrs.trna_added.renamed.gff3

#hapB
funannotate util gff-rename --gff3 Puccinia_striiformis_Pst104E_hapB.genes.utrs.trna_added.gff --fasta ../../../assembly_versions/v3.9_gapfill/v3.9.hapB.fasta --locus_tag Pst104E137 --out Puccinia_striiformis_Pst104E_hapB.genes.utrs.trna_added.renamed.gff3 --numbering 15639

# remove the Alias attribute in all lines
sed "s/Alias.*//g" Puccinia_striiformis_Pst104E_hapA.genes.utrs.trna_added.renamed.gff3 > Puccinia_striiformis_Pst104E_hapA.genes.utrs.trna_added.renamed.alias_rm.gff
sed "s/Alias.*//g" Puccinia_striiformis_Pst104E_hapB.genes.utrs.trna_added.renamed.gff3 > Puccinia_striiformis_Pst104E_hapB.genes.utrs.trna_added.renamed.alias_rm.gff3

