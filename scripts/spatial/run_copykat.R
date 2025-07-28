library(Seurat)
library(copykat)
library(ggplot2)
library(optparse)


option_list = list(
    make_option(c("-s", "--sample"), default=NULL, 
              help="sample name"),
    make_option(c("-c", "--cores"), type="integer", default=NULL, 
              help="cores", metavar="integer")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sample_name <- paste0(opt$sample, "_filtered")

dir.create(sample_name, showWarnings = FALSE)

# Prepare input matrix
exp.rawdata <- read.csv(paste0("../filtered_count_matrices/counts_", opt$sample, ".csv"), row.names=1)
rownames(exp.rawdata) <- sub("_.*", "", rownames(exp.rawdata))

# Define reference spots - spots annotated as non-tumor by pathologist or spots with the highest ESTIMATE score
healthy <- read.csv(paste0("chromothripsis/j462r/spatial_transcriptomics/CNVs/healthy_annotations/", opt$sample, "_healthy.csv"), row.names=1)
healthy <- healthy[row.names(healthy) %in% row.names(exp.rawdata), , drop=FALSE]
healthy_cells <- row.names(healthy)

exp.rawdata <- t(exp.rawdata)
print(colnames(exp.rawdata))
print(healthy_cells)

setwd(sample_name)

copykat_run <- copykat(
    rawmat=exp.rawdata, 
    id.type="S", 
    ngene.chr=5,
    # win size of 15 was used for the detection of ecDNA for samples B123 and B42
    win.size=25, 
    KS.cut=0.1, 
    sam.name=opt$sample, 
    distance="euclidean", 
    norm.cell.names=healthy_cells,
    output.seg="FALSE", 
    plot.genes="TRUE", 
    genome="hg20", 
    n.cores=opt$cores
)

pred <- data.frame(copykat_run$prediction)
pred <- pred[-which(pred$copykat.pred=="not.defined"),]
CNA <- data.frame(copykat_run$CNAmat)

saveRDS(pred, file = "predictions.Rda")
saveRDS(CNA, file = "CNA.Rda")
