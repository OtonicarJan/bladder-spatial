library(ComplexHeatmap)
library(RColorBrewer)
library(optparse)
library(circlize)
library(data.table)
library(viridis)


option_list = list(make_option(c("-s", "--sample"), default = NULL,
                               help = "sample name"))


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


sample_name = opt$sample

if (sample_name == "B154") {
    version = "_filtered_estimate/"
} else if (sample_name == "B178") {
    version = "_filtered_estimate2/" 
} else {
    version = "_filtered/"
}

copykat <-
  as.matrix(fread(
    paste0(
      "copykat/",
      sample_name,
      version,
      sample_name,
      "_copykat_bins.csv"
    )
  ), rownames = 1)
ont <-
  read.csv(paste0(
    "copykat/",
    sample_name,
    version,
    sample_name,
    "_ont_bins.csv"
  ))

predictions <-
  read.csv(
    paste0(
      "./copykat/",
      sample_name,
      version,
      sample_name,
      "_leiden_subclones_final.csv"
    ),
    row.names = 1
  )

segments <-
  read.csv(
    paste0(
      "./copykat/",
      sample_name,
      version,
      sample_name,
      "_spatialDE2_segments.csv"
    ),
    row.names = 1
  )

segments <- segments[!segments$segmentation_labels %in% c(8, 10, 11), , drop=FALSE]
aneuploid <- rownames(predictions[predictions$subclones != "diploid", , drop=FALSE])

aneuploid <- intersect(aneuploid, rownames(segments))
aneuploid_cells <- copykat[aneuploid,]

metadata <- read.csv(paste0("copykat/", sample_name, version, sample_name, "_metadata.csv"), row.names=1)
metadata <- metadata[rownames(aneuploid_cells), , drop = FALSE]
predictions <- segments[rownames(aneuploid_cells), , drop = FALSE]
print(predictions)

cnv_scores <- rowMeans(abs(aneuploid_cells))
cnv_scores_df <- data.frame(CNV_Score = cnv_scores)




col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
viridis_color_function <- colorRamp2(
  c(0, 0.05, 0.1),
  viridis(3)                  
)
viridis_color_function2 <- colorRamp2(
  c(-3000, 0, 3000),
  viridis(3)                  
)

cluster_colors <-
c(
  "0" = "#0173b2",
  "1" = "#de8f05",
  "2" = "#029e73",
  "3" = "#d55e00",
  "4" = "#cc78bc",
  "5" = "#ca9161",
  "6" = "#fbafe4",
    "7" = "#949494",
    "8" = "#56b4e9",
    "9" = "#ece133",
    "10" = "#99CBC3",
    "11" = "#CBC799"
)
row_ha = rowAnnotation(
    corr = metadata$correlation,
    cnv_score = cnv_scores_df$CNV_Score,
    estimate = metadata$estimate,
    col = list(
      corr = col_fun,
      cnv_score = viridis_color_function,
      estimate = viridis_color_function2
    ),
    annotation_name_side = "top",
    annotation_label = c("Correlation to bulk", "CNV score", "ESTIMATE score"),
    annotation_legend_param = list(
                            corr = list(title = "Correlation\nto bulk"),
                            cnv_score = list(title = "CNV score"),
                            estimate = list(title = "ESTIMATE score  ")
                        )
    )
row_cluster = rowAnnotation(
    cluster = predictions$segmentation_labels,
    col = list(
      cluster = cluster_colors
    ),
    annotation_legend_param = list(
                            cluster = list(title = "SpatialDE2\ncluster  ")
                        ),
    show_annotation_name = FALSE
)    

ont <- ont[ont$chr == "chr12", ]
ont <- ont[700:900, ]
s <- HeatmapAnnotation(ONT = anno_lines(ont[, ncol(ont)], axis_param = list(at = c(2, 4, 10), gp = gpar(fontsize = 10)), height = unit(5, "cm")), annotation_name_side = "left")

aneuploid_cells <- aneuploid_cells[ , grepl( "chr12" , colnames(aneuploid_cells))]
aneuploid_cells <- aneuploid_cells[ , 700:900]

col_fun2 = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))
h <- Heatmap(
    aneuploid_cells,
    name = "CN (log ratio)  ",
    col = col_fun2,
    cluster_columns = F,
    top_annotation = s,
    left_annotation = row_cluster,
    row_labels = F,
    use_raster = TRUE,
    row_split = predictions$segmentation_labels,
    row_gap = unit(0, "mm"),
    column_gap = unit(1, "mm"),
    border = T,
    show_row_names = F,
    show_row_dend = F,
    show_column_dend = F,
    show_column_names = F,
    row_title_rot = 0,
    column_title_rot = 45,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2"
)

svg(
    "fig4A_ecDNA_CNV.svg",
    width = 4,
    height = 6
)
draw(h)
dev.off()
