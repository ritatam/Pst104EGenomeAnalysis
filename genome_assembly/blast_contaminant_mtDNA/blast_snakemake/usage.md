
This script uses blastn locally. You will need a ncbi nt database downloaded locally on your machine.

first I'd recommend splitting your genome into smaller chunks, maximum 1Mbp per file. 

```python split_genome.py <input.fasta> <output_dir>```

set up configuration in the `config.yaml` file.

activate snakemake environment

```conda activate /home/groups/schwessinger/condaEnvs/snakemake```

install local blast environment

```snakemake --use-conda --conda-create-envs-only -c1```

dry-run to print commands to see if working properly

```snakemake --use-conda -np```

run the pipeline

```snakemake --use-conda -c8 --keep-going```

I included the flag "--keep-going" so the pipeline does not stop whenever it encounters problematic chunks that produce segmentation error "Segmentation fault (core dumped)". The reason for this error is unclear, though it might be due to nucleotide composition disliked by blast eg simple repeats (rarely). I did not include codes to automatically handle this error and had to further process manually. AFTER the pipeline finishes running, user must go to log folder in .snakemake and check for problematic chunks using `grep Error <logfile> -A3 | grep ".fasta"`. For these I simply picked them out and split them further into even smaller chunks (e.g. 30kbp).