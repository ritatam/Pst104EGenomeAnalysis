ref=ref/v3.9.hapA.fasta
query=ref/v3.9.hapB.fasta

l=200 
b=500
c=500

nucmer --maxmatch -l $l -b $b -c $c $ref $query 

delta-filter -m -i 90 -l 100 out.delta > out_m_i90_l100.delta
show-coords -THrd out_m_i90_l100.delta > out_m_i90_l100.coords

mkdir -p syri_i90_l100
syri -c out_m_i90_l100.coords -r $ref -q $query -d out_m_i90_l100.delta --dir syri_i90_l100

echo -e  "$ref\thapA\tlc:firebrick;lw:2\n$query\thapB\tlc:royalblue;lw:2" > genomes.txt
rm chrorder.txt; for i in {1..18}; do echo -e "chr${i}A" >> chrorder.txt; done

plotsr --sr syri_i90_l100/syri.out --genomes genomes.txt -o syri_i90_l100_plot.png --chrord chrorder.txt -H 5 -W 6 -d 600 -f 5
