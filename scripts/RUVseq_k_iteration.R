# Run DESeq2 for all comparisons on exacloud
library("dplyr")
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("ggrepel")

args <- commandArgs()

help <- function(){
  cat("run DESeq2 for TCGA and GTEx data\n")
  cat("Usage: \n")
  cat("--dds     : dds object                                                         [ required ]\n")
  cat("--res   : results file                                                            [ required ]\n")
  cat("--comparison   : Comparison (target-vs-baseline)                                             [ default = 0.25 ]\n")
   cat("\n")
  q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
  help()
} else {
  dds        <- sub('--dds=', '', args[grep('--dds=', args)])
  res        <- sub( '--res=', '',args[grep('--res=',args)])
  contrast   <- sub( '--contrast=', '',args[grep('--contrast=',args)])
  kvalue     <- sub( '--kvalue=', '',args[grep('--kvalue=',args)])
}

target <- unlist(strsplit(contrast, "-"))[1]
baseline_TCGA <- unlist(strsplit(contrast, "-"))[3]
baseline_GTEx <- unlist(strsplit(contrast, "-"))[4]

baseline <- gsub("_GTEx","",baseline_GTEx)
target <- gsub("_TCGA","",target)

outDir <- paste0("results/", contrast, "_RUVseq")

upCol = "#FF9999"
downCol = "#99CCFF"
ncCol = "#CCCCCC"
adjp <- 0.01
FC <- 2

# Read in dds object and results file
dds <- readRDS(dds)
res <- read.table(res, sep = "\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

res <- res[order(res$padj, decreasing = TRUE),]

kvalue <- as.numeric(kvalue)

# Take top 5000 least DE genes
genes <- rownames(res)[1:5000]

# Generate EDAseq object
library(EDASeq)
library(RColorBrewer)
library(RUVSeq)
df <- as.data.frame(colData(dds))
colors <- c("firebrick1","black")

library(plyr)
df$Color <- plyr::mapvalues(df$Group, unique(df$Group), colors)

set <- newSeqExpressionSet(as.matrix(counts(dds)), phenoData = df)

pdf(paste0(outDir, "/", contrast, "-boxplot-initial.pdf"))
plotRLE(set, outline = FALSE, ylim = c(-4,4), col=as.vector(df$Color))
dev.off()

pdf(paste0(outDir, "/", contrast, "-PCA-initial.pdf"))
plotPCA(set, col=as.vector(df$Color))
dev.off()

set <- betweenLaneNormalization(set, which = "upper")

pdf(paste0(outDir, "/", contrast, "-boxplot-initial-UQ.pdf"))
plotRLE(set, outline = FALSE, ylim = c(-4,4), col=as.vector(df$Color))
dev.off()

pdf(paste0(outDir, "/", contrast, "-PCA-initial-UQ.pdf"))
plotPCA(set, col=as.vector(df$Color))
dev.off()

setk <- RUVg(x=set, cIdx=genes, k=kvalue)

pdf(paste0(outDir, "/", contrast, "-boxplot-UQ-k-", kvalue, ".pdf"))
plotRLE(setk, outline = FALSE, ylim = c(-4,4), col=as.vector(df$Color))
dev.off()

pdf(paste0(outDir, "/", contrast, "-PCA-UQ-k-", kvalue, ".pdf"))
plotPCA(setk, col=as.vector(df$Color))
dev.off()

dds <- DESeqDataSetFromMatrix(countData = counts(setk), colData = pData(setk), design = ~ W_1 + Type)

dds <- DESeq(dds)
res <- results(dds, contrast=c("Type", target, baseline),  independentFiltering = FALSE, cooksCutoff = Inf)
res <- res[order(res$padj),]

write.table(as.data.frame(res), paste0(outDir, "/", contrast, "-results-RUVseqNorm-k", kvalue, ".txt"))

saveRDS(dds, paste0(outDir, "/", contrast, "-dds-RUVseqNorm-k", kvalue, ".rds"

