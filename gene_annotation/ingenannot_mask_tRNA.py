import pandas as pd
import os

gff="Puccinia_striiformis_Pst104E_hapA.gff3"
trna_masked_gff="Puccinia_striiformis_Pst104E_hapA.trna_commented.gff3"
df = pd.read_csv(gff,sep="\t",skiprows=1, header=None)

trna_list = [n.split(";")[0][3:][:-3] for n in df[df[2] == "tRNA"][8].tolist()]
with open(gff, "r") as file:
    lines = file.readlines()
newlines = []
for line in lines:
    comment_line = any(trna in line for trna in trna_list)
    newlines.append("#" + line.rstrip() if comment_line else line.rstrip())
with open(trna_masked_gff, "w") as out:
    for n in newlines:
        print(n, file=out)
