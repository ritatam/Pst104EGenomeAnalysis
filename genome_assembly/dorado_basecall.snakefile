# snakemake pipeline for duplex basecalling with dorado on National Computational Infrastructure, ANU, Australia.
# see https://github.com/ritatam/dorado-snake-nci/tree/main.

configfile: "config.yaml"

from glob import glob
OUTPUT_DIR = config["output_dir"]
POD5 = glob(f"{config['pod5']}/*.pod5")

rule all:
    input:
        f"{OUTPUT_DIR}/pairs_from_bam/",
        f"{OUTPUT_DIR}/pod5_splitduplex/",
        f"{OUTPUT_DIR}/split_duplex_pair_ids.txt",
        f"{OUTPUT_DIR}/final/duplex_orig.fastq",
        f"{OUTPUT_DIR}/final/duplex_splitduplex.fastq",
        f"{OUTPUT_DIR}/final/split_on_adapters"
        


rule simplex_basecall:
    input:
        list(POD5)
    output:
        f"{OUTPUT_DIR}/unmapped_reads_with_moves.sam"
    envmodules:
        "cuda/11.7.0"
    shell:
        '{config[dorado]} basecaller --device "cuda:0,1" {config[model]} {config[pod5]} --verbose --emit-moves > {output}'

rule samtools_sam_to_bam:
    input:
        f"{OUTPUT_DIR}/unmapped_reads_with_moves.sam"
    output:
        bam = f"{OUTPUT_DIR}/unmapped_reads_with_moves.bam",
        bai = f"{OUTPUT_DIR}/unmapped_reads_with_moves.bam.bai"
    envmodules:
        "samtools/1.12"
    shell:
        "samtools view -@8 -S -b {input} > {output.bam} && "
        "samtools index -@8 {output.bam}"

rule find_duplex_pairs:
    input:
        bam = f"{OUTPUT_DIR}/unmapped_reads_with_moves.bam",
        bai = f"{OUTPUT_DIR}/unmapped_reads_with_moves.bam.bai"
    output:
        pairs_dir = directory(f"{OUTPUT_DIR}/pairs_from_bam/"),
        pair_ids_filtered = f"{OUTPUT_DIR}/pairs_from_bam/pair_ids_filtered.txt"
    params:
        duplex_tools = config["duplex_tools_env"]
    shell:
        ". {params.duplex_tools} && "
        "duplex_tools pair --output_dir {output.pairs_dir} {input.bam}"

rule find_duplex_pairs_from_nonsplit_reads:
    input:
        sam = f"{OUTPUT_DIR}/unmapped_reads_with_moves.sam",
        pod5 = config["pod5"]
    output:
        pod5_splitduplex = directory(f"{OUTPUT_DIR}/pod5_splitduplex/"),
        split_pair_txt = f"{OUTPUT_DIR}/split_duplex_pair_ids.txt"
    params:
        duplex_tools = config["duplex_tools_env"]
    shell:
        ". {params.duplex_tools} && "
        "duplex_tools split_pairs --threads 24 {input.sam} {input.pod5} {output.pod5_splitduplex} && "
        "cat {output.pod5_splitduplex}/*_pair_ids.txt > {output.split_pair_txt}"

rule stereo_basecall_from_main_pairing:
    input:
        f"{OUTPUT_DIR}/pairs_from_bam/pair_ids_filtered.txt"
    output:
        f"{OUTPUT_DIR}/final/duplex_orig.fastq"
    envmodules:
        "cuda/11.7.0"
    shell:
        '{config[dorado]} duplex {config[model]} {config[pod5]} --device "cuda:0,1" --threads 48 --pairs {input} --emit-fastq > {output}'

rule stereo_basecall_from_additional_pairing:
    input:
        f"{OUTPUT_DIR}/split_duplex_pair_ids.txt"
    output:
        f"{OUTPUT_DIR}/final/duplex_splitduplex.fastq"
    params:
        f"{OUTPUT_DIR}/pod5_splitduplex/"
    envmodules:
        "cuda/11.7.0"
    shell:
        '{config[dorado]} duplex {config[model]} {params} --threads 48 --device "cuda:0,1" --pairs {input} --emit-fastq > {output}'

rule split_chimeric_reads_on_adapter:
    input:
        f"{OUTPUT_DIR}/final/duplex_orig.fastq",
        f"{OUTPUT_DIR}/final/duplex_splitduplex.fastq"
    output:
        directory(f"{OUTPUT_DIR}/final/split_on_adapters"),
    params:
        duplex_tools = config["duplex_tools_env"],
        fastq_dir = f"{OUTPUT_DIR}/final"
    shell:
        ". {params.duplex_tools} && "
        "duplex_tools split_on_adapter --threads 24 {params.fastq_dir} {output} Native"