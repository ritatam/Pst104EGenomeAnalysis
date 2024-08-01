library(bambu)
library(doParallel)
library(foreach)
library(dplyr)
 

fa <- 'v3.9.fasta'
# using CDS here to avoid overlapping exon features from spurious gene annotation.
# in this file, "CDS" is replaced by "exon" in the third column
# grep -P "\tCDS\t Puccinia_striiformis_Pst104E.gtf | sed -i "s/CDS/exon/g" Puccinia_striiformis_Pst104E.CDS.gtf
gtf <- 'Puccinia_striiformis_Pst104E.CDS.gtf'
anno <- prepareAnnotations(gtf)
bam_list <- c('Pst104E_U_UG_rep1.ont.mapped.bam', 'Pst104E_U_UG_rep2.ont.mapped.bam',
              'Pst104E_U_UG_rep3.ont.mapped.bam', 'Pst104E_U_UG_rep4.ont.mapped.bam',
              'Pst104E_U_4dpi_rep1.ont.mapped.bam', 'Pst104E_U_4dpi_rep2.ont.mapped.bam',
              'Pst104E_U_4dpi_rep3.ont.mapped.bam', 'Pst104E_U_4dpi_rep4.ont.mapped.bam',
              'Pst104E_U_6dpi_rep1.ont.mapped.bam', 'Pst104E_U_6dpi_rep2.ont.mapped.bam',
              'Pst104E_U_6dpi_rep3.ont.mapped.bam', 'Pst104E_U_6dpi_rep4.ont.mapped.bam',
              'Pst104E_U_8dpi_rep1.ont.mapped.bam', 'Pst104E_U_8dpi_rep2.ont.mapped.bam',
              'Pst104E_U_8dpi_rep3.ont.mapped.bam', 'Pst104E_U_8dpi_rep4.ont.mapped.bam',
              'Pst104E_U_10dpi_rep1.ont.mapped.bam', 'Pst104E_U_10dpi_rep2.ont.mapped.bam',
              'Pst104E_U_10dpi_rep3.ont.mapped.bam', 'Pst104E_U_10dpi_rep4.ont.mapped.bam',
              'Pst104E_U_12dpi_rep1.ont.mapped.bam', 'Pst104E_U_12dpi_rep2.ont.mapped.bam',
              'Pst104E_U_12dpi_rep3.ont.mapped.bam', 'Pst104E_U_12dpi_rep4.ont.mapped.bam')
bam_list <- paste("bam/", bam_list, sep = "")

bambu_quant <- function(bam) {
  sample_name <- gsub(".ont.mapped.bam", "", basename(bam))
  se <- bambu(reads=bam, annotations=anno, genome=fa, discovery=FALSE, opt.discovery = list(min.txScore.singleExon = 0))  # Use 1 core for each bambu call
  writeBambuOutput(se, path = file.path("./bambu_quant", sample_name))
}

# parallelise
num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

foreach(bam=bam_list, .packages = c("bambu")) %dopar% {
  bambu_quant(bam)
}

# stop parallel backend
stopCluster(cl)

# now join the output count_genes.txt dataframes together.
merged_counts <- data.frame()
for (bam in bam_list) {
    sample_name <- gsub(".ont.mapped.bam", "", basename(bam))
    col_name <- gsub(".bam", "", basename(bam))
    counts_gene_txt <- file.path("./bambu_quant", sample_name, "counts_gene.txt")
    counts_data <- read.table(counts_gene_txt, header=TRUE, stringsAsFactors=F)
    counts_subset <- counts_data[, c("GENEID", col_name)]
    if (nrow(merged_counts) == 0) {
        merged_counts <- counts_subset
    } else {
        merged_counts <- left_join(merged_counts, counts_subset, by = "GENEID")
    }
}

output_path <- "./all_cond.gene_counts.cds.bambu.tsv"
write.table(merged_counts, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("merged counts have been written to", output_path, "\n")
