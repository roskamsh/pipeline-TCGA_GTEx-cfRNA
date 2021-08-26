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
  cat("--dds     : counts file                                                         [ required ]\n")
  cat("--rld   : metadata                                                            [ required ]\n")
  cat("--contrast   : Comparison (target-vs-baseline)                                             [ default = 0.25 ]\n")
   cat("\n")
  q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args)) ){
  help()
} else {
  dds     <- sub('--dds=', '', args[grep('--dds=', args)])
  rld   <- sub( '--rld=', '',args[grep('--rld=',args)])
  contrast   <- sub( '--contrast=', '',args[grep('--contrast=',args)])
}

target <- unlist(strsplit(contrast, "-"))[1]
baseline <- unlist(strsplit(contrast, "-"))[3]

outDir <- paste0("results/", contrast)

dds <- readRDS(dds)
rld <- readRDS(rld)

upCol = "#FF9999"
downCol = "#99CCFF"
ncCol = "#CCCCCC"
adjp <- 0.01
FC <- 2

# Pairwise PCA Plot
pdf(paste0(outDir, "/", contrast, "-PCA.pdf"))
plotPCA(rld, intgroup="Type")
dev.off()

res <- results(dds, contrast=c("Type", target, baseline),  independentFiltering = FALSE, cooksCutoff = Inf)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=c("Type", target, baseline), res=res)

# MA plot - calc norm values yourself to plot with ggplot
# MA plot is log2normalized counts (averaged across all samples) vs. log2FC

# extract normalized counts to calculate values for MA plot
norm_counts <- counts(dds, normalized=TRUE)

## select up regulated genes
forPlot <- as.data.frame(res)
forPlot$log2Norm <- log2(rowMeans(norm_counts))
forPlot$Gene <- rownames(forPlot)

up <- forPlot$padj < adjp & forPlot$log2FoldChange > log2(FC)
sum(up)

## select down regulated genes
down <- forPlot$padj < adjp & forPlot$log2FoldChange < -log2(FC)
sum(down)

# Grab the top 5 up and down regulated genes to label in the volcano plot
if (sum(up)>5) {
  temp <- forPlot[up,]
  upGenesToLabel <- head(rownames(temp[order(-temp$log2FoldChange),], 5))
} else if (sum(up) %in% 1:5) {
  temp <- forPlot[up,]
  upGenesToLabel <- rownames(temp[order(-temp$log2FoldChange),])
}

if (sum(down)>5) {
  temp <- forPlot[down,]
  downGenesToLabel <- head(rownames(temp[order(temp$log2FoldChange),], 5))
} else if (sum(down) %in% 1:5) {
  temp <- forPlot[down,]
  downGenesToLabel <- rownames(temp[order(temp$log2FoldChange),])
}

forPlot$Expression <- ifelse(down, 'down',
                             ifelse(up, 'up','NS'))
forPlot$Expression <- factor(forPlot$Expression, levels=c("up","down","NS"))

# Assign colours to conditions
if (sum(up)==0 & sum(down)==0) {
  colours <- ncCol
} else if (sum(up)==0) {
  colours <- c(downCol, ncCol)
} else if (sum(down)==0) {
  colours <- c(upCol, ncCol)
} else {
  colours <- c(upCol, downCol, ncCol)
}

# Create vector for labelling the genes based on whether genes are DE or not
if (exists("downGenesToLabel") & exists("upGenesToLabel")) {
  genesToLabel <- c(downGenesToLabel, upGenesToLabel)
} else if (exists("downGenesToLabel") & !exists("upGenesToLabel")) {
  genesToLabel <- downGenesToLabel
} else if (!exists("downGenesToLabel") & exists("upGenesToLabel")) {
  genesToLabel <- upGenesToLabel
}

if (exists("genesToLabel")) {
  maPlot <- ggplot(forPlot, mapping=aes(x=log2Norm, y=log2FoldChange, colour=Expression)) +
    geom_point() +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="black") +
    geom_label_repel(aes(label=ifelse(Gene %in% genesToLabel, as.character(Gene),'')),box.padding=0.1, point.padding=0.5, segment.color="gray70", show.legend=FALSE) +
    scale_colour_manual(values=colours) +
    ggtitle(paste(baseline, "vs", target)) +
    xlab("log2(Normalized counts)") +
    ylab("log2(Fold Change)") +
    theme(plot.title = element_text(hjust=0.5))
} else {
  maPlot <- ggplot(forPlot, mapping=aes(x=log2Norm, y=log2FoldChange, colour=Expression)) +
    geom_point() +
    geom_hline(yintercept=c(-1,1), linetype="dashed", color="black") +
    scale_colour_manual(values=colours) +
    ggtitle(paste(baseline, "vs", target)) +
    xlab("log2(Normalized counts)") +
    ylab("log2(Fold Change)") +
    theme(plot.title = element_text(hjust=0.5))
}

# MA plot
pdf(paste0(outDir, "/", contrast, "-ma-plot.pdf"))
print({
  maPlot
})
dev.off()

res <- res[order(res$padj),]
write.table(as.data.frame(res), file=paste0(outDir, "/", contrast, "-DEG.txt"), quote=FALSE, sep='\t')

# Heatmap of top 50 genes
topGenes <- head(order(res$padj),50)

df <- as.data.frame(colData(rld))

annot <- as.data.frame(cbind(rownames(df), paste(df$Type)))
names(annot) <- c("SampleID", "Type")
rownames(annot) <- annot$SampleID
annot$SampleID <- NULL

hm <- pheatmap(assay(rld), cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of all genes:", baseline, "vs", target), color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

save_pheatmap_pdf <- function(x, filename) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(hm, paste0(outDir, "/", contrast, "-heatmap-all-genes.pdf"))

hm <- pheatmap(assay(rld)[topGenes,], cluster_rows=T, scale="row", fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T, cluster_cols=T, annotation_col=annot, labels_col=as.character(rownames(df)), main = paste("Heatmap of top 50 DE genes:", baseline, "vs", target), color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

save_pheatmap_pdf(hm, paste0(outDir, "/", contrast, "-heatmap-top50-DE-genes.pdf"))

