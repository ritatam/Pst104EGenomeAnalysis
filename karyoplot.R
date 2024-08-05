library(biomaRt)
library(regioneR)
library(karyoploteR)
library(tidyr)
library(dplyr)
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE)

# genome
Genome <- read.csv("v3.9.chr.haplotype-paired.fasta.fai", sep="\t", header=FALSE)[,1:2]
Genome$start <- 0
colnames(Genome) <- c("chr", "end", "start")
Genome <- Genome[, c("chr", "start", "end")]
# centromere
centromeres <- read.csv("Pst104E_centromeres_refined.csv")
# genes
Genes <- read.csv("v3.9.genedensity", sep="\t", header=FALSE)[,c(1:3,7)]
Genes[,2] <- Genes[,2]+1 # GRanges uses 1-based coordinates.
colnames(Genes) <- c("chr", "start", "end", "genedensity")
# TEs
TEs <- read.csv("TEs/v3.9_REPET_TEs.repeatdensity", sep="\t", header=FALSE)[,c(1:3,7)]
colnames(TEs) <- c("chr", "start", "end", "TEdensity")
# assembly gaps
gaps <- read.csv("Pst104E_assemblyGaps.csv", header=FALSE)
colnames(gaps) <- c("chr", "start", "end")
# rDNA arrays
rDNA_arrays <- data.frame(chr=c("chr13A","chr13B"), start=c(1915422,1855615), end=c(1994741,1914026))
# secretome and effectors
seceff <- read.csv("Puccinia_striiformis_Pst104E.secretome-effectors.TM_anno_cleaned.gff3", sep="\t", header=FALSE)[,c(1,4,5,8)]
colnames(seceff) <- c("chr", "start", "end", "name")
# mating type genes from Yan
matingtype <- read.csv("Pst104E_MatingTypeAnno_Yan_edited.gff3", sep="\t", header=FALSE)[,c(1,4,5,9)]
colnames(matingtype) <- c("chr", "start", "end", "name")


# plot parameters settings
pp <- getDefaultPlotParams(plot.type=1)
pp$data1height <- 180
pp$data1outmargin <- 150
pp$ideogramheight <- 50
pp$leftmargin <- 0.1
pp$bottommargin <- 800


plot <- function(haplotype) {
  
  # load all data 
  cpgmeth <- read.table(paste0("v3.9.10k.windows.Pst104E.duplex.5mCG-5hmCG.hap", haplotype, ".segmeth.tsv"), sep="\t", header=TRUE)[,2:10]
  colnames(cpgmeth) <- c("chr", "start", "end", "window_name", "strand",
                         "m_meth_calls", "m_unmeth_calls", "m_no_calls", "m_methfrac")
  Genome_hap <- toGRanges(subset(Genome, substr(chr, nchar(chr), nchar(chr)) == haplotype))
  Genes_hap <- toGRanges(subset(Genes, substr(chr, nchar(chr), nchar(chr)) == haplotype))
  TEs_hap <- subset(TEs, substr(chr, nchar(chr), nchar(chr)) == haplotype)
  cpgmeth_hap <- subset(cpgmeth, substr(chr, nchar(chr), nchar(chr)) == haplotype)
  centromeres_hap <- subset(centromeres, substr(chr, nchar(chr), nchar(chr)) == haplotype)
  matingtype_hap <- subset(matingtype, substr(chr, nchar(chr), nchar(chr)) == haplotype)
  
  # lay out the chromosomes
  kp <- plotKaryotype(genome=Genome_hap, plot.type=1, labels.plotter=NULL, plot.params=pp)
  # add chromosome names
  kpAddChromosomeNames(kp, xoffset=0.02, yoffset=0.5,cex=1,adj=c(0.9,NA), font=2)
  # add base numbers and ticks
  kpAddBaseNumbers(kp, tick.dist = 500000, tick.len=40, units="Mb", add.units="Mb", cex=1,
                   minor.tick.dist=100000, minor.tick.len=25, minor.tick.col="#00000050")
  
  ###################
  ### data tracks ###
  ###################
  # Track 1: gene density
  if (haplotype=="A") {
    genedensity_cols <- c("#FFF4F4", "#EFADB5", "#780000")}
  if (haplotype=="B") {
    genedensity_cols <- c("#F4F5FF", "#ADB1EF", "#030A72")}
  kp <- kpHeatmap(kp, data=Genes_hap, y=Genes_hap$genedensity, col=genedensity_cols, data.panel = "ideogram")
  # Track 2: TE density
  at <- autotrack(current.track=c(2,3), total.tracks=6)
  chr_list_hap <- unique(TEs_hap$chr)
  kpArea(kp, chr=TEs_hap$chr, x=(TEs_hap$start+TEs_hap$end)/2, y=TEs_hap$TEdensity, r0=at$r0, r1=at$r1, col="#A992AD",border=NA)
  # Track 3: methylation density
  at <- autotrack(current.track=c(4,5), total.tracks=6)
  for (chr in chr_list_hap) {
    chr_data <- cpgmeth_hap[cpgmeth_hap$chr == chr, ]
    # kpArea does weird stuff with polygons so a quick dirty fix is to add a 1bp 
    # window at min and max coordinates and assign methfrac value as zero such 
    # that polygon base is fixed at y=0
    min_coord <- min(chr_data$start)
    max_coord <- max(chr_data$end)
    min_row <- data.frame(chr = chr,start = min_coord,end = min_coord + 1,window_name = "NoName",strand = ".",m_meth_calls = 0,m_unmeth_calls = 0,m_no_calls = 0,m_methfrac = 0)
    max_row <- data.frame(chr = chr,
                          start = max_coord,end = max_coord + 1,window_name = "NoName",strand = ".",m_meth_calls = 0,m_unmeth_calls = 0,m_no_calls = 0,m_methfrac = 0)
    chr_data <- rbind(min_row, chr_data, max_row)
    chr_data <- chr_data[order(chr_data$start), ]
    chr_data <- replace_na(chr_data, replace = list(m_methfrac = 0))
    kpPolygon(kp, chr=chr_data$chr, x=(chr_data$start+chr_data$end)/2, y=chr_data$m_methfrac, r0=at$r0, r1=at$r1, col="#333A73", border=NA)
  }
  
  ######################
  ### other features ###
  ######################
  at <- autotrack(current.track=1, total.tracks=6)
  kpRect(kp, chr=rDNA_arrays$chr, x0=rDNA_arrays$start, x1=rDNA_arrays$end, y0=-0.4, y1=0.95, r0=at$r0, r1=at$r1, col="darkgreen", border=NA)
  # gaps
  kpPoints(kp, chr=gaps$chr, x=(gaps$end-gaps$start)/2+gaps$start, y=0.05, data.panel=1, bg="red", cex=0.8, pch=21, lwd=1, col="#FFFFFF")
  # chromosome boundary
  kpRect(kp, chr=Genome$chr, x0=Genome$start, x1=Genome$end, y0=0, y1=1, data.panel="ideogram")
  # secretome/effector genes
  kpRect(kp, chr=seceff$chr, x0=seceff$start, x1=seceff$end, y0=-0.2, y1=0.8, r0=at$r0, r1=at$r1, col="#FF7A00", border=NA)
  # centromeres 
  kpRect(kp, chr=centromeres_hap$chr, x0=centromeres_hap$start, x1=centromeres_hap$end, y0=0, y1=0.9, data.panel="all", border=NA, col="#00000027")
  # mating type genes
  kpPlotMarkers(kp, chr=matingtype_hap$chr, x=(matingtype_hap$end-matingtype_hap$start)/2+matingtype_hap$start, labels=matingtype_hap$name, text.orientation="horizontal", label.dist=0.005, y=0.9, cex=0.82, font=4)
  # telomere 
  for (chr in unique(Genome$chr)) {
    chr_data <- Genome[Genome$chr == chr, ]
    # p
    kpPoints(kp, chr=chr_data$chr, x=0, y=0.2, r0=at$r0, r1=at$r1, col="black", pch=25, cex=1.2, bg="black")
    # q
    if (chr != "chr16A"){
      kpPoints(kp, chr=chr_data$chr, x=chr_data$end, y=0.2, r0=at$r0, r1=at$r1, col="black", pch=25, cex=1.2, bg="black")
    }
  }
  
  ####################
  ###### legend ######
  ####################
  legend(x=0.7, y=0.3, inset=.02, pt.bg=c("#333A73", "#A992AD", "#00000000", "#00000050", "#FF7A00","black", "darkgreen", "red"), 
         legend=c("CpG methylation density", "TE density", "", "Centromere", "Secretome gene","Telomere",  "Ribosomal DNA array", "Gap"),
         pch = c(22, 22, 22,    22,    22, 25, 22, 21),
         col=c(NA, NA, NA,   NA, NA, NA, NA, NA),
         pt.cex=c(2,2,0,   2,2,1.3,2,1.5))
}


# plot karyoplot
png(filename="haplotypeA_karyoplot_v2.png",height=20*300, width=10*300, res=300)
kp_A <- plot("A")
dev.off()

png(filename="haplotypeB_karyoplot_v2.png",height=20*300, width=10*300, res=300)
kp_B <- plot("B")
dev.off()


# plot gene density colour bar separately
png(filename="gene_density_scale_hapA.png",height=2.5*300, width=2*300, res=300)
genedensity_cols <- c("#FFF4F4", "#780000")
values <- seq(0, 1, length.out = 100)
par(mar = c(5, 3, 2, 6)) # bottom, left, top, right margins
z <- matrix(values, nrow = 1)
image(1, values, z, col = colorRampPalette(genedensity_cols)(100), 
      axes = FALSE, xlab = "", ylab = "")
axis(4, at = c(0, 1), labels = c("0",  "1"), las = 3, padj=-0.6)
mtext("Gene density", side = 2, line = 1, padj=1)
dev.off()

png(filename="gene_density_scale_hapB.png",height=2.5*300, width=2*300, res=300)
genedensity_cols <- c("#F4F5FF", "#030A72")
values <- seq(0, 1, length.out = 100)
par(mar = c(5, 3, 2, 6)) # bottom, left, top, right margins
z <- matrix(values, nrow = 1)
image(1, values, z, col = colorRampPalette(genedensity_cols)(100), 
      axes = FALSE, xlab = "", ylab = "")
axis(4, at = c(0, 1), labels = c("0",  "1"), las = 3, padj=-0.6)
mtext("Gene density", side = 2, line = 1, padj=1)
dev.off()
