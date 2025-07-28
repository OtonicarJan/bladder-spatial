library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)


txtFontSize=12
axisFontSize=12
axisTtlFontSize=16
lgdTtlFontSize=16
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), 
                   legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), 
                   strip.background = element_blank(), axis.line=element_line(linewidth=0.7, color="black"), 
                   axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), 
                   axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), 
                   legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), 
                   text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))

cov_ILL <-
  read.table(
    "../../../../bladder_cancer/liquid_biopsy/bam_files/B123/ichorCNA_B123.cna.seg",
    header = TRUE,
    sep = "\t"
  )
colnames(cov_ILL) <-
  c("chr",
    "start",
    "end",
    "copy_number",
    "event",
    "logR",
    "subclone_status",
    "corrected_copynumber",
    "corrected_call",
    "CN")

cov_ILL$logR <- cov_ILL$logR - median(cov_ILL$logR, na.rm = TRUE)
cov_ILL <- cov_ILL[!is.na(cov_ILL$CN), ]

cov_ILL[cov_ILL$CN > 28, ] <- 28
cov_ILL$start <- cov_ILL$start - 1


cov_ONT <-
  read.table(
    "../../../../bladder_cancer/ont/cnvs/coverage/B123tumor.cov.gz",
    header = TRUE,
    sep = "\t"
  )

colnames(cov_ONT) <-
  c("chr",
    "start",
    "end",
    "mappable",
    "counts",
    "CN")

event_colors <- c(
  "GAIN" = "#FF8263",
  "AMP" = "#FF643D",
  "HLAMP" = "#FF4518",
  "HLAMP2" = "#FF3200",
  "HETD" = "#377eb8",
  "NEUT" = "gray50"
)


plot_list <- list()
widths <- numeric()

for (chrlabel in unique(cov_ONT$chr)) {
  
  ont_chr <- cov_ONT[cov_ONT$chr == chrlabel, ]
  ill_chr <- cov_ILL[cov_ILL$chr == chrlabel, ]
  
  # midpoints
  mid_ont <- (ont_chr$start + ont_chr$end) / 2
  mid_ill <- (ill_chr$start + ill_chr$end) / 2
  
  df_ont <- data.frame(pos = mid_ont, signal = ont_chr$CN)
  df_ill <- data.frame(pos = mid_ill, signal = ill_chr$logR, event = ill_chr$corrected_call)
  
  chr_width <- max(c(df_ont$pos, df_ill$pos), na.rm = TRUE) - min(c(df_ont$pos, df_ill$pos), na.rm = TRUE)
  widths <- c(widths, chr_width)
  
  p_ont <- ggplot(df_ont, aes(x = pos, y = signal)) +
    geom_point(size = 0.2) + ylab(NULL) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "gray30", linewidth = 0.3) +
    scale_x_continuous(labels=NULL) +
    scale_y_continuous(labels=NULL, breaks = seq(0, 28, by = 4), limits = c(0, 28)) +
    scienceTheme +
    theme(
      axis.title.x = element_blank()
    )
  
  p_ill <- ggplot(df_ill, aes(x = pos, y = signal, color = event)) +
    geom_point(size = 1.2) + ylab(NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.3) +
    scale_x_continuous(
      labels=NULL,
      limits = c(min(c(df_ont$pos, df_ill$pos), na.rm = TRUE), max(c(df_ont$pos, df_ill$pos), na.rm = TRUE))
    ) +
    scale_color_manual(values = event_colors) +
    scale_y_continuous(labels=NULL, breaks = seq(-0.5, 1, by = 0.25), limits = c(-0.5, 1)) +
    xlab(chrlabel) +
    scienceTheme +
    theme(legend.position = "none")
  
  p_combined <- plot_grid(p_ont, p_ill, ncol = 1, align = "v", rel_heights = c(1, 1))
  plot_list[[chrlabel]] <- p_combined
  
}

rel_widths <- widths / sum(widths)
final_plot <- plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "h", rel_widths = rel_widths)

ggsave(file="Fig2AB_B123_ONT_cfDNA.png", width=40, height=12, limitsize = FALSE)



# Zoom on chromosome 6
plot_list <- list()
widths <- numeric()

for (chrlabel in c("chr6")) {
  
  ont_chr <- cov_ONT[cov_ONT$chr == chrlabel, ]
  ill_chr <- cov_ILL[cov_ILL$chr == chrlabel, ]
  
  # midpoints
  mid_ont <- (ont_chr$start + ont_chr$end) / 2
  mid_ill <- (ill_chr$start + ill_chr$end) / 2
  
  df_ont <- data.frame(pos = mid_ont, signal = ont_chr$CN)
  df_ill <- data.frame(pos = mid_ill, signal = ill_chr$logR, event = ill_chr$corrected_call)
  
  chr_width <- max(c(df_ont$pos, df_ill$pos), na.rm = TRUE) - min(c(df_ont$pos, df_ill$pos), na.rm = TRUE)
  widths <- c(widths, chr_width)
  
  p_ont <- ggplot(df_ont, aes(x = pos, y = signal)) +
    geom_point(size = 0.2) + ylab(NULL) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "gray30", linewidth = 0.3) +
    scale_x_continuous(labels=NULL) +
    scale_y_continuous(labels=NULL, breaks = seq(0, 28, by = 4), limits = c(0, 28)) +
    scienceTheme +
    theme(
      axis.title.x = element_blank()
    )
  
  p_ill <- ggplot(df_ill, aes(x = pos, y = signal, color = event)) +
    geom_point(size = 0.5) + ylab(NULL) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.3) +
    scale_x_continuous(
      labels=NULL,
      limits = c(min(c(df_ont$pos, df_ill$pos), na.rm = TRUE), max(c(df_ont$pos, df_ill$pos), na.rm = TRUE))
    ) +
    scale_color_manual(values = event_colors) +
    scale_y_continuous(labels=NULL, breaks = seq(-0.5, 1, by = 0.25), limits = c(-0.5, 1)) +
    xlab(chrlabel) +
    scienceTheme +
    theme(legend.position = "none")
  
  p_combined <- plot_grid(p_ont, p_ill, ncol = 1, align = "v", rel_heights = c(1, 1))
  plot_list[[chrlabel]] <- p_combined
  
}



# Final chromosome panel (horizontally arranged)
rel_widths <- widths / sum(widths)
final_plot <- plot_grid(plotlist = plot_list, ncol = length(plot_list), align = "h", rel_widths = rel_widths)

ggsave(file="Fig2AB_B123_chr6_zoom.png", width=4, height=4)
