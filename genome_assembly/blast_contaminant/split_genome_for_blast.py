# splits input genome assembly fasta file into multiple fasta files each containing 1Mbp for subsequent blastn.
# usage: python split_genome.py <input.fasta> <output_dir>

import pandas as pd
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math


def main():
    parser = argparse.ArgumentParser(description="Split an input fasta file into multiple fasta files, each containing 1Mbp.")
    parser.add_argument("input_fasta", help="Path to the input fasta file.")
    parser.add_argument("output_dir", help="Path to the output directory where the split fasta files will be saved.")
    args = parser.parse_args()

    input_fasta = args.input_fasta
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    split_fasta_per_1mbp(input_fasta, output_dir)


def split_fasta_per_1mbp(fasta_fn, outdir):
    '''
    Split an input fasta (e.g. genome assembly with many chromosome-scale contigs) into many fasta files, each containing 1Mbp.
    For example, if a sequence in the input fasta is 2843686bp, it will be split into three fasta files,
    named <sequence_id>_{n}.fasta. The first two will contain 1,000,000bp, and the last one will contain 843686bp.
    '''
    if os.path.exists(outdir) == False:
        os.makedir(outdir)
    outpaths = []
    limit = 1000000
    for record in SeqIO.parse(fasta_fn, "fasta"):
        seq_len = len(record.seq)
        print(f"Length of {record.id}: {seq_len}bp")
        if seq_len > limit:
            num_chunks = math.ceil(seq_len/limit)
            print(f"Splitting {record.id} into {num_chunks} chunks.")
            for i in range(1,num_chunks+1):
                chunk = str(record.seq[(i-1)*limit:i*limit])
                print(f"Chunk {i-1}: {len(chunk)}bp.")
                chunk_fn = os.path.join(outdir, f"{record.id}_{i-1}.fasta")
                newrecord = SeqRecord(Seq(chunk), id=f"{record.id}_{i-1}", description="")
                SeqIO.write(newrecord, chunk_fn, "fasta")
                outpaths.append(chunk_fn)
        else:
            chunk_fn = os.path.join(outdir, f"{record.id}.fasta")
            SeqIO.write(record, chunk_fn, "fasta")
            outpaths.append(chunk_fn)
        print(" ")
    return outpaths

if __name__ == "__main__":
    main()