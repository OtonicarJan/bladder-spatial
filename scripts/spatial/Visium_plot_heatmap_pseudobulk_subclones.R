library(ComplexHeatmap)
library(RColorBrewer)
library(optparse)
library(circlize)
library(data.table)
library(viridis)

# First plot for the sample B123
copykat <- as.data.frame(fread("copykat/subclones_pseudobulks.csv"), rownames=1)

copykat <- copykat[copykat$sample == "B123", ]
copykat <- copykat[copykat$subclone != "diploid", ]

metadata <- as.data.frame(copykat[ ,c("subclone", "sample")])
metadata$sample <- factor(metadata$sample, levels = samples)

chromosomes <- sub(":.*", "", colnames(copykat))
chromosomes <- factor(chromosomes, levels = paste0("chr", c(1:22, 'X')))

col_fun2 = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))
h <- Heatmap(
    as.matrix(copykat),
    col = col_fun2,
    name = "CN (log ratio)  ",
    column_split = chromosomes,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    cluster_rows = TRUE,
    row_labels = FALSE,
    use_raster = TRUE,
    column_gap = unit(0, "mm"),
    row_gap = unit(0, "mm"),
    border = TRUE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    column_title_rot = 90,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
  )

svg(paste0("fig2E_subclones_heatmap_B123.svg"), width=8, height=2)
draw(h)
dev.off()



# Plot for all samples
samples <- c("B22", "B24", "B60", "B154", "B156", "B175", "B178", "B4", "B42", "B123")
copykat <- as.data.frame(fread("copykat/subclones_pseudobulks.csv"), rownames=1)

copykat <- copykat[copykat$subclone != "diploid", ]

print(colnames(copykat))

metadata <- as.data.frame(copykat[ ,c("subclone", "sample")])
metadata$sample <- factor(metadata$sample, levels = samples)

copykat <- copykat[ , !(names(copykat) %in% c("V1", "subclone", "sample"))]


chromosomes <- sub(":.*", "", colnames(copykat))
chromosomes <- factor(chromosomes, levels = paste0("chr", c(1:22, 'X')))

str(chromosomes)

col_fun2 = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))
h <- Heatmap(
    as.matrix(copykat),
    col = col_fun2,
    name = "CN (log ratio)  ",
    column_split = chromosomes,
    cluster_columns = FALSE,
    row_split = metadata$sample,
    cluster_row_slices = FALSE,
    cluster_rows = TRUE,
    row_labels = FALSE,
    use_raster = TRUE,
    column_gap = unit(0, "mm"),
    row_gap = unit(0, "mm"),
    border = TRUE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    column_title_rot = 45,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    show_parent_dend_line = FALSE,
    row_title_gp = gpar(fontsize = 16),
    row_title_rot = 0
  )

svg(paste0("suppfig_6A_subclones_heatmap.svg"), width=12, height=8)
draw(h)
dev.off()
