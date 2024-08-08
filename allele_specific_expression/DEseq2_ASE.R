
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(readr)
library(pasilla)
library(vsn)
library(pheatmap)

create_dir <- function(dirname) {
  if (file.exists(dirname)) {
    cat("This directory already exists\n")
  } else {
    dir.create(dirname)
  }
}

DEanalysis_results <- function(dds, query, ref) {
  DEres_name <- paste("condition", query, "vs", ref, sep="_")
  DEres <- results(dds, name=DEres_name)
  DEres <- results(dds, contrast=c("condition", query, ref))
  DEres <- DEres[order(DEres$padj), ] # order rows by padj
  write.csv(as.data.frame(DEres), file=paste(paste(DEres_csv_outpath, DEres_name, sep="/"), ".csv", sep=""))
  plotMA(DEres, ylim=c(-10, 10), main=paste(DEres_name, "MA-plot results", sep=" "))
  return(DEres)
}

csv2df <- function(comparison) {
  fn <- paste(DEres_csv_outpath, comparison, sep="/")
  fn <- paste(fn, ".tsv", sep="")
  df <- read_csv(file=fn, col_types=cols())
  return(df)
}


create_dir("./DESeq_results_table")
create_dir("./DESeq_LFCshrinkage_results_table")
DEres_csv_outpath <- "./DESeq_results_table"
DELFCres_csv_outpath <- "./DESeq_LFCshrinkage_results_table"


# read in bambu count matrix already processed by bambu_quant_processing script to meet deseq2 format requirement
counts <- as.matrix(read.csv("allele_paired.gene_counts.cds.FILTERED.bambu.tsv", sep="\t", row.names="allele_pair_id"))
# make sure columns of the count matrix and the rows of column data are in the same order
coldata <- data.frame(condition = factor(rep(c("UG_a", "4dpi_a", "6dpi_a", "8dpi_a", "10dpi_a", "12dpi_a", "UG_b", "4dpi_b", "6dpi_b", "8dpi_b", "10dpi_b", "12dpi_b"), each=4)))
rownames(coldata) <- colnames(counts)
coldata

# count matrix already pre-filtered for abundance>=5 in bambu_quant_processing script
# now convert count matrix to DESeqDataSet format
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)


# run DEseq2 w/ hapA genes as reference
Run_DESeq2 <- function(dds, query, ref) {
  dds$condition <- relevel(dds$condition, ref=ref)
  dds <- DESeq(dds)
  results <- DEanalysis_results(dds, query, ref)
  results <- results[!is.na(results$padj) & results$padj < 0.05, ]
  print(paste("number of results for", query,"_vs_",ref,": ", nrow(results)))
  return(results)
}
UG <- Run_DESeq2(dds, "UG_b", "UG_a")
dpi4 <- Run_DESeq2(dds, "4dpi_b", "4dpi_a")
dpi6 <- Run_DESeq2(dds, "6dpi_b", "6dpi_a")
dpi8 <- Run_DESeq2(dds, "8dpi_b", "8dpi_a")
dpi10 <- Run_DESeq2(dds, "10dpi_b", "10dpi_a")
dpi12 <- Run_DESeq2(dds, "12dpi_b", "12dpi_a")


# save normalised data matrix for other use. note DESeq2 doesn't use normalised counts.
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds <- estimateSizeFactors(dds)
# normalised counts transformation
normalised_counts <- counts(dds, normalized=TRUE)
write.table(normalised_counts, file="allele_paired.bambu_DESeq2_norm_counts.csv", quote=F)


