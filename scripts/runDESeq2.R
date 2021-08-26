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
  cat("--counts     : counts file                                                         [ required ]\n")
  cat("--metadata   : metadata                                                            [ required ]\n")
  cat("--comparison   : Comparison (target-vs-baseline)                                             [ default = 0.25 ]\n")
   cat("\n")
  q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
  help()
} else {
  counts        <- sub('--counts=', '', args[grep('--counts=', args)])
  metadata        <- sub( '--metadata=', '',args[grep('--metadata=',args)])
  contrast   <- sub( '--contrast=', '',args[grep('--contrast=',args)])
}

target <- unlist(strsplit(contrast, "-"))[1]
baseline <- unlist(strsplit(contrast, "-"))[3]

outDir <- paste0("results/", contrast)

dir.create(outDir, recursive = TRUE)

# Read in metadata table and order according to sampleID
md <- read.delim(file=metadata, sep = "\t", stringsAsFactors = FALSE)
md <- md[order(md$SampleID),]

# Read in counts table
subdata <- read.table(counts, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
subdata <- subdata[,order(colnames(subdata))]

# Extract only the Types that we want in further analysis & only the PP_ID and Status informative columns
md <- filter(md, Type == baseline | Type == target)

# Keep only the PP_IDs of the types we have chosen in the metadata table above
rownames(md) <- md$SampleID
md$SampleID <- NULL
keep <- colnames(subdata)[colnames(subdata) %in% rownames(md)]
subdata <- subdata[, keep]
dim(subdata)

# Check
stopifnot(rownames(md)==colnames(subdata))

# Round data
subdata$gene_id <- NULL
rounded <- round(subdata, 0)

# Obtain the number of genes that meet padj<0.01 for reference line in histogram
dds <- DESeqDataSetFromMatrix(countData=rounded,
                              colData=md,
                              design= ~ Type)

dds <- estimateSizeFactors(dds)

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) >= 1, ]

# Normalization and pre-processing
dds <- DESeq(dds)

saveRDS(dds, file=paste0(outDir, "/", contrast, "_dds.rds"))

# obtain normalized counts
rld <- vst(dds, blind=FALSE)
saveRDS(rld, file=paste0(outDir, "/", contrast, "_rlog.rds"))

