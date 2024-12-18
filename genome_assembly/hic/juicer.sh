# 1. prepare reference
mkdir references
cp assembly.v3.9.fasta references/
bwa index assembly.v3.9.fasta

# 2. prepare fastq
mkdir fastq
cp hic_fastq.gz fastq

# 3. softlink core juicer scripts 
ln -s /media/nvme/rita/softwares/juicer/CPU scripts

# 4. generate restriction site positions
wget https://github.com/aidenlab/juicer/blob/main/misc/generate_site_positions.py
# generate_site_positions.py edits:
# patterns = {'DpnII-HinFI-MseI-DdeI': ['GATC', 'GANTC', 'TTAA', 'CTNAG']}
# filenames = {'Pst104E_v3.9': '/media/nvme/rita/project/104e_verkko_v2/juicer/references/assembly.v3.9.fasta'}
mkdir restriction_sites
python generate_site_positions.py DpnII-HinFI-MseI-DdeI Pst104E_v3.9

# 5. prepare chromosome size file
bioawk -c fastx '{print $name"\t"length($seq)}' assembly.v3.9.fasta > ../chrom.sizes

# 6. used custom python script to generate ligation sites; see generate_ligation_site.py.
# then edited juicer.sh
# Set ligation junction based on restriction enzyme
# if [ -z "$ligation" ]
# then
# ...
#        DpnII-HinFI-MseI-DdeI) ligation="'(GATCGATC|GAATAATC|GAATATTC|GAATACTC|GAATAGTC|GATTAATC|GATTATTC|GATTACTC|GATTAGTC|GACTAATC|GACTATTC|GACTACTC|GACTAGTC|GAGTAATC|GAGTATTC|GAGTACTC|GAGTAGTC|TTATAA|CTAATAAG|CTAATTAG|CTAATCAG|CTAATGAG|CTTATAAG|CTTATTAG|CTTATCAG|CTTATGAG|CTCATAAG|CTCATTAG|CTCATCAG|CTCATGAG|CTGATAAG|CTGATTAG|CTGATCAG|CTGATGAG)'" ;;      <------
#        none) ligation="XXXX";;
# ...

# 7. run juicer
export wdir=/media/nvme/rita/project/104e_verkko_v2/juicer
bash ./scripts/juicer.sh -S early -D ${wdir} -d ${wdir} -g Pst104E_v3.9 -z ${wdir}/references/assembly.v3.9.fasta -y ${wdir}/restriction_sites/Pst104E_v2_DpnII-HinFI-MseI-DdeI.txt -p ${wdir}/chrom.sizes -t 8 -b DpnII-HinFI-MseI-DdeI