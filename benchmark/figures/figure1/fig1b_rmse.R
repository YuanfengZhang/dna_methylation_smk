# 2025-0820
# EQA unified_rmse
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

rmse <- read.table("unified_rmse.csv",
  sep = ",", header = T
)
head(rmse)
rmse1 <- rmse %>%
  bind_rows(
    rmse %>%
      filter(lab == "NP1") %>%
      filter(grepl("_1$", sample)) %>% # 找出所有 *_1 的行
      mutate(sample = sub("_1$", "_2", sample)) # 把 _1 改成 _2
  ) %>%
  select(lab, sample, rmse)
head(rmse1)

rmse1_data <- rmse1 %>%
  pivot_wider(names_from = sample, values_from = rmse)
head(rmse1_data)
rmse1_data <- as.data.frame(rmse1_data)
rownames(rmse1_data) <- rmse1_data$lab
rmse1_data <- rmse1_data[, -1]
# anno_20
anno_20 <- as.data.frame(rownames(rmse1_data))
colnames(anno_20) <- "Filename"
anno_20$Platform <- sapply(strsplit(as.character(anno_20$Filename), ""), function(x) {
  paste(x[1], x[2], sep = "")
})

anno_20$Platform <- gsub("BS", "WGBS", anno_20$Platform)
anno_20$Platform <- gsub("EM", "EM-seq", anno_20$Platform)
anno_20$Platform <- gsub("GM", "GM-seq", anno_20$Platform)
anno_20$Platform <- gsub("RR", "RRBS", anno_20$Platform)
anno_20$Platform <- gsub("RM", "RRMS", anno_20$Platform)
anno_20$Platform <- gsub("NP", "Nanopore", anno_20$Platform)

rownames(anno_20) <- anno_20$Filename
# anno_20<-anno_20[,-1]
anno_20$Filename <- NULL
# anno_12
anno_12 <- as.data.frame(colnames(rmse1_data))
colnames(anno_12) <- "Filename"
anno_12$Type <- sapply(strsplit(as.character(anno_12$Filename), ""), function(x) {
  paste(x[1], x[2], sep = "")
})
rownames(anno_12) <- anno_12$Filename
# anno_12<-anno_12[,-1]
anno_12$Filename <- NULL

ann_colors3 <- list(
  Platform = c(
    WGBS = "#FFDDAD",
    "EM-seq" = "#AEC5EB",
    RRBS = "#049A8F",
    "GM-seq" = "#3A405A",
    RRMS = "#BEE5A0",
    Nanopore = "#A983C6",
    MA = "#3C5487"
  ),
  Type = c(
    D5 = "#4CC3D9",
    D6 = "#7BC8A4",
    F7 = "#FFC65D",
    M8 = "#F16745",
    T1 = "gray100",
    T2 = "gray66",
    T3 = "gray33",
    T4 = "gray0",
    BC = "#69001f",
    BL = "#f7b293"
  )
)

lab_order <- c(
  "MA2", "MA1", "MA3", "GM3", "RR1", "GM1", "BS1", "GM2", "BS3",
  "EM2", "EM4", "EM3", "RM1", "BS2", "EM1", "NP1",
  "BS4"
)
sample_order <- c(
  "D5_1", "D5_2", "D6_1", "D6_2", "F7_1", "F7_2", "M8_1", "M8_2",
  "T1_1", "T1_2", "T2_1", "T2_2", "T3_1", "T3_2", "T4_1", "T4_2",
  "BC_1", "BC_2", "BL_1", "BL_2"
)
rmse1_data1 <- rmse1_data[lab_order, sample_order]
rmse_p1 <- pheatmap(rmse1_data1,
  cluster_cols = F,
  cluster_rows = F,
  annotation_legend = T,
  treeheight_row = 16, treeheight_col = 16,
  scale = "none",
  # scale = "column",
  # gaps_row = c(12),
  gaps_col = c(8, 16),
  breaks = seq(9, 27, by = 2),
  color = viridis(9, option = "viridis"),
  cellwidth = 14, cellheight = 14,
  clustering_distance_cols = "correlation",
  clustering_distance_rows = "correlation",
  # clustering_distance_cols = "euclidean",
  # clustering_distance_rows = "euclidean",
  show_colnames = F, show_rownames = T,
  annotation_row = anno_20,
  annotation_col = anno_12,
  #  border=F,
  border_color = "black",
  annotation_colors = ann_colors3
)
rmse_p1

topptx(rmse_p1, "rmse_heatmap.pptx", width = 8, height = 6)

max(rmse1_data1)
min(rmse1_data1)

# protocol_sample
rmse_p <- rmse1 %>%
  mutate(protocol = substr(lab, 1, 2))
rmse_p$protocol <- gsub("BS", "WGBS", rmse_p$protocol)
rmse_p$protocol <- gsub("EM", "EM-seq", rmse_p$protocol)
rmse_p$protocol <- gsub("GM", "GM-seq", rmse_p$protocol)
rmse_p$protocol <- gsub("RR", "RRBS", rmse_p$protocol)
rmse_p$protocol <- gsub("RM", "RRMS", rmse_p$protocol)
rmse_p$protocol <- gsub("NP", "Nanopore", rmse_p$protocol)

rmse_p$protocol <- factor(rmse_p$protocol,
  levels = c("WGBS", "EM-seq", "RRBS", "GM-seq", "RRMS", "Nanopore", "MA")
)
rmse_p$lab <- factor(rmse_p$lab,
  levels = c(
    "BS1", "BS2", "BS3", "BS4",
    "EM1", "EM2", "EM3", "EM4",
    "RR1",
    "GM1", "GM2", "GM3",
    "RM1",
    "NP1",
    "MA1", "MA2", "MA3"
  )
)
rmse_p$sample <- factor(rmse_p$sample,
  levels = c("D5", "D6", "F7", "M8", "T1", "T2", "T3", "T4", "BC", "BL")
)
method_colors <- c(
  "D5" = "#4CC3D9", "D6" = "#7BC8A4", "F7" = "#FFC65D", "M8" = "#F16745",
  "BC" = "#69001f", "BL" = "#f7b293",
  "T1" = "gray100", "T2" = "gray66", "T3" = "gray33", "T4" = "gray0",
  "HF" = "#8039f9",
  "BS1" = "#FFDDAD", "BS2" = "#FFDDAD", "BS3" = "#FFDDAD", "BS4" = "#FFDDAD",
  "EM1" = "#AEC5EB", "EM2" = "#AEC5EB", "EM3" = "#AEC5EB", "EM4" = "#AEC5EB",
  "RR1" = "#049A8F",
  "GM1" = "#3A405A", "GM2" = "#3A405A", "GM3" = "#3A405A",
  "RM1" = "#BEE5A0",
  "NP1" = "#A983C6",
  "MA1" = "#3C5487", "MA2" = "#3C5487", "MA3" = "#3C5487"
)
protocol_colors <- c(
  "WGBS" = "#FFDDAD",
  "EM-seq" = "#AEC5EB",
  "RRBS" = "#049A8F",
  "GM-seq" = "#3A405A",
  "RRMS" = "#BEE5A0",
  "Nanopore" = "#A983C6",
  "MA" = "#3C5487"
)
p <- ggplot(
  data = rmse_p,
  aes(x = lab, y = rmse, fill = lab)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(
    side = "r", errorbar.draw = FALSE,
    outlier.size = 1,
    outlier.stroke = 0.2,
    width = 0.2, linewidth = 0.5
  ) +
  geom_half_point_panel(aes(fill = sample),
    side = "l",
    range_scale = .85,
    shape = 21, size = 1.5, color = "black"
  ) +
  scale_fill_manual(values = method_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black")
  ) +
  facet_grid2(. ~ protocol,
    scales = "free", space = "free_x",
    strip = strip_nested(
      background_x = elem_list_rect(fill = protocol_colors, color = NA),
      text_x = element_text(color = "black", face = "bold")
    )
  ) +
  labs(y = "RMSE", x = "Lab")
p
topptx(p, "protocol_rmse.pptx", width = 8, height = 4)
