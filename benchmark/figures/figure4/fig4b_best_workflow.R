# 2025-0822
# EQA best_workflow
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(tidyverse)
library(pheatmap)
library(readxl)
library(grid)
library(gridExtra)
library(patchwork)
library(cowplot)

file_path <- "best_workflow.xlsx"

aligner_file <-  read_excel(file_path, sheet = 1)
counter_file <-  read_excel(file_path, sheet = 2)
deduper_file <-  read_excel(file_path, sheet = 3)
calibrator_file <-  read_excel(file_path, sheet = 4)
# BS_EM
# aligner
plot_BS_aligner <- function(aligner_file, metrics, breaks, color){
  aligner <- aligner_file %>%
    filter(method == "BS_EM") %>%
    select(tools_aligner, metrics) %>%
    as.data.frame()
  rownames(aligner) <- aligner$tools_aligner
  aligner <- aligner[-1]
  
  p <- pheatmap(aligner,
                cluster_cols = F,
                cluster_rows = F,
                annotation_legend = T,
                treeheight_row = 6,treeheight_col = 7,
                scale = "none",
                # scale = "column",
                # gaps_row = c(12),
                # gaps_col = c(24, 42),
                breaks = breaks,
                color = color,
                cellwidth = 14, cellheight = 14,
                clustering_distance_cols="correlation",
                clustering_distance_rows="correlation",
                #clustering_distance_cols = "euclidean",
                #clustering_distance_rows = "euclidean",
                show_colnames =F,show_rownames = T,
                #annotation_row=anno_20,
                #annotation_col=anno_counter,
                #  border=F,
                border_color = "black")
  g <- grid::grid.grabExpr(print(p))
  return(g)
}

p_D6 <- plot_BS_aligner(aligner_file, "D6_rmse", breaks = seq(16.5, 18.5, by = 0.4), color = colorRampPalette(c("#FFC65D","white",'gray33'))(6))
p_BC <- plot_BS_aligner(aligner_file, "BC_rmse", breaks = seq(15, 17.5, by = 0.5), color = colorRampPalette(c("#FFC65D","white",'gray33'))(6))
p_aligner <- plot_grid(p_D6, p_BC, ncol = 1, align = "v")
p_aligner

topptx(p_aligner,"BS_aligner.pptx",width = 11,height = 7)

align <-  aligner_file %>%
  filter(method == "BS_EM")
min(align$D6_rmse, na.rm = TRUE)
max(align$D6_rmse, na.rm = TRUE)
min(align$BC_rmse, na.rm = TRUE)
max(align$BC_rmse, na.rm = TRUE)

# counter
plot_BS_counter <- function(counter_file, metrics, breaks, color){
  counter <- counter_file %>%
    filter(method == "BS_EM") %>%
    select(tools_aligner, tools_counter, metrics) %>%
    pivot_wider(names_from = tools_aligner, values_from = metrics) %>%
    as.data.frame()
  rownames(counter) <- counter$tools_counter
  counter <- counter[-1]
  
  anno_counter <- data.frame(aligner = colnames(counter))
  rownames(anno_counter) <- colnames(counter)
  
  ann_colors3= list(
    aligner=c('bwa-meth'='#E64B3599',bsmapz="#4DBBD599",'batmeth2'="#00A08799",'bismark-hisat2'="#3C548899",'bismark-bowtie2'="#F39B7F99",fame="#8491B499"))  
  
  num_labels <- counter
  num_labels[] <- ifelse(is.na(counter), "NA", "")
  
  p <- pheatmap(counter,
                         cluster_cols = F,
                         cluster_rows = F,
                         annotation_legend = T,
                         treeheight_row = 6,treeheight_col = 7,
                         scale = "none",
                         # scale = "column",
                         # gaps_row = c(12),
                         # gaps_col = c(24, 42),
                         breaks = breaks,
                         color = color,
                         cellwidth = 14, cellheight = 14,
                         clustering_distance_cols="correlation",
                         clustering_distance_rows="correlation",
                         #clustering_distance_cols = "euclidean",
                         #clustering_distance_rows = "euclidean",
                         show_colnames =F,show_rownames = T,
                         #annotation_row=anno_20,
                         annotation_col=anno_counter,
                         #  border=F,
                         border_color = "black",annotation_colors = ann_colors3,
                         display_numbers = num_labels,   
                         number_color = "white",
                         na_col = "#3C548899")
  g <- grid::grid.grabExpr(print(p))
  return(g)
}

p_D6 <- plot_BS_counter(counter_file, "D6_rmse", breaks = seq(16, 36, by = 2), color = colorRampPalette(c("#FFC65D","white",'gray33'))(11))
p_BC <- plot_BS_counter(counter_file, "BC_rmse", breaks = seq(15, 39, by = 2), color = colorRampPalette(c("#FFC65D","white",'gray33'))(13))
p_counter <- plot_grid(p_D6, p_BC, ncol = 1, align = "v")
p_counter

topptx(p_counter,"BS_counter.pptx",width = 11,height = 7)

count <- counter_file %>%
  filter(method == "BS_EM")
min(count$D6_rmse, na.rm = TRUE)
max(count$D6_rmse, na.rm = TRUE)
min(count$BC_rmse, na.rm = TRUE)
max(count$BC_rmse, na.rm = TRUE)
# deduper
plot_BS_deduper <- function(deduper_file, metrics, breaks, color){
  deduper <- deduper_file %>%
    filter(method == "BS_EM") %>%
    mutate(tools_aligner_counter = paste(tools_aligner, tools_counter, sep = "_")) %>%
    select(tools_aligner_counter, tools_deduper, metrics) %>%
    pivot_wider(names_from = tools_aligner_counter, values_from = metrics) %>%
    as.data.frame()
  rownames(deduper) <- deduper$tools_deduper
  deduper <- deduper[-1]
  
  anno_deduper <- data.frame(aligner_counter = colnames(deduper))
  rownames(anno_deduper) <- colnames(deduper)
  
  ann_colors3= list(
    aligner_counter=c('bwa-meth_astair'="#91D1C299",'bwa-meth_methyldackel'="#DC000099",'bsmapz_astair'="#7E614899",'bsmapz_methyldackel'="#B09C8599"))  
     
  num_labels <- deduper
  num_labels[] <- ifelse(is.na(deduper), "NA", "")
  
  p <- pheatmap(deduper,
                cluster_cols = F,
                cluster_rows = F,
                annotation_legend = T,
                treeheight_row = 6,treeheight_col = 7,
                scale = "none",
                # scale = "column",
                # gaps_row = c(12),
                # gaps_col = c(24, 42),
                breaks = breaks,
                color = color,
                cellwidth = 14, cellheight = 14,
                clustering_distance_cols="correlation",
                clustering_distance_rows="correlation",
                #clustering_distance_cols = "euclidean",
                #clustering_distance_rows = "euclidean",
                show_colnames =F,show_rownames = T,
                #annotation_row=anno_20,
                annotation_col=anno_deduper,
                #  border=F,
                border_color = "black",
                annotation_colors = ann_colors3,
                display_numbers = num_labels,   
                number_color = "white",
                na_col = "#3C548899")
  g <- grid::grid.grabExpr(print(p))
  return(g)
}

p_D6 <- plot_BS_deduper(deduper_file, "D6_rmse", breaks = seq(16.5, 17.5, by = 0.1), color = colorRampPalette(c("#FFC65D","white",'gray33'))(11))
p_BC <- plot_BS_deduper(deduper_file, "BC_rmse", breaks = seq(15.2, 15.6, by = 0.1), color = colorRampPalette(c("#FFC65D","white",'gray33'))(5))
p_deduper <- plot_grid(p_D6, p_BC, ncol = 1, align = "v")
p_deduper

topptx(p_deduper,"BS_deduper.pptx",width = 11,height = 7)

dedup <- deduper_file %>%
  filter(method == "BS_EM")
min(dedup$D6_rmse, na.rm = TRUE)
max(dedup$D6_rmse, na.rm = TRUE)
min(dedup$BC_rmse, na.rm = TRUE)
max(dedup$BC_rmse, na.rm = TRUE)

# calibrator
plot_BS_calibrator <- function(calibrator_file, metrics, breaks, color){
  calibrator <- calibrator_file %>%
    filter(method == "BS_EM") %>%
    mutate(tools_aligner_counter_dedup_cali = paste(tools_aligner, tools_counter, tools_deduper, tools_calibrator, sep = "_")) %>%
    select(tools_aligner_counter_dedup_cali, metrics) %>%
    as.data.frame()
  rownames(calibrator) <- calibrator$tools_aligner_counter_dedup_cali
  calibrator <- calibrator[-1]
  
  p <- pheatmap(calibrator,
                cluster_cols = F,
                cluster_rows = F,
                annotation_legend = T,
                treeheight_row = 6,treeheight_col = 7,
                scale = "none",
                # scale = "column",
                # gaps_row = c(12),
                # gaps_col = c(24, 42),
                breaks = breaks,
                color = color,
                cellwidth = 14, cellheight = 14,
                clustering_distance_cols="correlation",
                clustering_distance_rows="correlation",
                #clustering_distance_cols = "euclidean",
                #clustering_distance_rows = "euclidean",
                show_colnames =F,show_rownames = T,
                #annotation_row=anno_20,
                #annotation_col=anno_counter,
                #  border=F,
                border_color = "black")
  g <- grid::grid.grabExpr(print(p))
  return(g)
}

p_D6 <- plot_BS_calibrator(calibrator_file, "D6_rmse", breaks = seq(16.5, 35, by = 2), color = colorRampPalette(c("#FFC65D","white",'gray33'))(9))
p_BC <- plot_BS_calibrator(calibrator_file, "BC_rmse", breaks = seq(15, 23, by = 1), color = colorRampPalette(c("#FFC65D","white",'gray33'))(9))
p_calibrator <- plot_grid(p_D6, p_BC, ncol = 1, align = "v")
p_calibrator

topptx(p_calibrator,"BS_calibrator.pptx",width = 11,height = 7)

cali <-  calibrator_file %>%
  filter(method == "BS_EM")
min(cali$D6_rmse, na.rm = TRUE)
max(cali$D6_rmse, na.rm = TRUE)
min(cali$BC_rmse, na.rm = TRUE)
max(cali$BC_rmse, na.rm = TRUE)


"#E64B354C" "#4DBBD54C" "#00A0874C" "#3C54884C" "#F39B7F4C" "#8491B44C" "#91D1C24C" "#DC00004C" "#7E61484C" "#B09C854C"
"#E64B3599" "#4DBBD599" "#00A08799" "#3C548899" "#F39B7F99" "#8491B499" "#91D1C299" "#DC000099" "#7E614899" "#B09C8599"
pal_npg("nrc", alpha=0.6)(10)
pal_npg("nrc", alpha=0.9)(10)
"#E64B35E5" "#4DBBD5E5" "#00A087E5" "#3C5488E5" "#F39B7FE5" "#8491B4E5" "#91D1C2E5" "#DC0000E5" "#7E6148E5" "#B09C85E5"

# GM
# aligner
plot_GM_aligner <- function(aligner_file, metrics, breaks, color){
  aligner <- aligner_file %>%
    filter(method == "GM") %>%
    select(tools_aligner, metrics) %>%
    as.data.frame()
  rownames(aligner) <- aligner$tools_aligner
  aligner <- aligner[-1]
  
  p <- pheatmap(aligner,
                cluster_cols = F,
                cluster_rows = F,
                annotation_legend = T,
                treeheight_row = 6,treeheight_col = 7,
                scale = "none",
                # scale = "column",
                # gaps_row = c(12),
                # gaps_col = c(24, 42),
                breaks = breaks,
                color = color,
                cellwidth = 14, cellheight = 14,
                clustering_distance_cols="correlation",
                clustering_distance_rows="correlation",
                #clustering_distance_cols = "euclidean",
                #clustering_distance_rows = "euclidean",
                show_colnames =F,show_rownames = T,
                #annotation_row=anno_20,
                #annotation_col=anno_counter,
                #  border=F,
                border_color = "black")
  g <- grid::grid.grabExpr(print(p))
  return(g)
}

p_D6 <- plot_GM_aligner(aligner_file, "D6_rmse", breaks = seq(14, 15, by = 0.1), color = colorRampPalette(c("#FFC65D","white",'gray33'))(11))
p_BC <- plot_GM_aligner(aligner_file, "BC_rmse", breaks = seq(14, 15.5, by = 0.3), color = colorRampPalette(c("#FFC65D","white",'gray33'))(6))
p_aligner <- plot_grid(p_D6, p_BC, ncol = 1, align = "v")
p_aligner

topptx(p_aligner,"GM_aligner.pptx",width = 11,height = 7)

align <-  aligner_file %>%
  filter(method == "GM")
min(align$D6_rmse, na.rm = TRUE)
max(align$D6_rmse, na.rm = TRUE)
min(align$BC_rmse, na.rm = TRUE)
max(align$BC_rmse, na.rm = TRUE)

# counter
plot_GM_counter <- function(counter_file, metrics, breaks, color){
  counter <- counter_file %>%
    filter(method == "GM") %>%
    select(tools_aligner, tools_counter, metrics) %>%
    pivot_wider(names_from = tools_aligner, values_from = metrics) %>%
    as.data.frame()
  rownames(counter) <- counter$tools_counter
  counter <- counter[-1]
  
  anno_counter <- data.frame(aligner = colnames(counter))
  rownames(anno_counter) <- colnames(counter)
  
  ann_colors3= list(
    aligner=c('gem3'='#E64B3599',bowtie2="#4DBBD599",'biscuit'="#00A08799"))  
  
  num_labels <- counter
  num_labels[] <- ifelse(is.na(counter), "NA", "")
  
  p <- pheatmap(counter,
                cluster_cols = F,
                cluster_rows = F,
                annotation_legend = T,
                treeheight_row = 6,treeheight_col = 7,
                scale = "none",
                # scale = "column",
                # gaps_row = c(12),
                # gaps_col = c(24, 42),
                breaks = breaks,
                color = color,
                cellwidth = 14, cellheight = 14,
                clustering_distance_cols="correlation",
                clustering_distance_rows="correlation",
                #clustering_distance_cols = "euclidean",
                #clustering_distance_rows = "euclidean",
                show_colnames =F,show_rownames = T,
                #annotation_row=anno_20,
                annotation_col=anno_counter,
                #  border=F,
                border_color = "black",annotation_colors = ann_colors3,
                display_numbers = num_labels,   
                number_color = "white",
                na_col = "#3C548899")
  g <- grid::grid.grabExpr(print(p))
  return(g)
}

p_D6 <- plot_GM_counter(counter_file, "D6_rmse", breaks = seq(13.5, 27.5, by = 2), color = colorRampPalette(c("#FFC65D","white",'gray33'))(8))
p_BC <- plot_GM_counter(counter_file, "BC_rmse", breaks = seq(14, 33, by = 2), color = colorRampPalette(c("#FFC65D","white",'gray33'))(11))
p_counter <- plot_grid(p_D6, p_BC, ncol = 1, align = "v")
p_counter

topptx(p_counter,"GM_counter.pptx",width = 11,height = 7)

count <- counter_file %>%
  filter(method == "GM")
min(count$D6_rmse, na.rm = TRUE)
max(count$D6_rmse, na.rm = TRUE)
min(count$BC_rmse, na.rm = TRUE)
max(count$BC_rmse, na.rm = TRUE)

# deduper
plot_GM_deduper <- function(deduper_file, metrics, breaks, color){
  deduper <- deduper_file %>%
    filter(method == "GM") %>%
    mutate(tools_aligner_counter = paste(tools_aligner, tools_counter, sep = "_")) %>%
    select(tools_aligner_counter, tools_deduper, metrics) %>%
    pivot_wider(names_from = tools_aligner_counter, values_from = metrics) %>%
    as.data.frame()
  rownames(deduper) <- deduper$tools_deduper
  deduper <- deduper[-1]
  
  anno_deduper <- data.frame(aligner_counter = colnames(deduper))
  rownames(anno_deduper) <- colnames(deduper)
  
  ann_colors3= list(
    aligner_counter=c('gem3_rastair'="#3C548899",'gem3_methyldackel'="#F39B7F99",'bowtie2_methyldackel'="#8491B499",'bowtie2_rastair'="#91D1C299"))  
     
  num_labels <- deduper
  num_labels[] <- ifelse(is.na(deduper), "NA", "")
  
  p <- pheatmap(deduper,
                cluster_cols = F,
                cluster_rows = F,
                annotation_legend = T,
                treeheight_row = 6,treeheight_col = 7,
                scale = "none",
                # scale = "column",
                # gaps_row = c(12),
                # gaps_col = c(24, 42),
                breaks = breaks,
                color = color,
                cellwidth = 14, cellheight = 14,
                clustering_distance_cols="correlation",
                clustering_distance_rows="correlation",
                #clustering_distance_cols = "euclidean",
                #clustering_distance_rows = "euclidean",
                show_colnames =F,show_rownames = T,
                #annotation_row=anno_20,
                annotation_col=anno_deduper,
                #  border=F,
                border_color = "black",
                annotation_colors = ann_colors3,
                display_numbers = num_labels,   
                number_color = "white",
                na_col = "#3C548899")
  g <- grid::grid.grabExpr(print(p))
  return(g)
}

p_D6 <- plot_GM_deduper(deduper_file, "D6_rmse", breaks = seq(13.5, 14.5, by = 0.1), color = colorRampPalette(c("#FFC65D","white",'gray33'))(11))
p_BC <- plot_GM_deduper(deduper_file, "BC_rmse", breaks = seq(14, 14.5, by = 0.1), color = colorRampPalette(c("#FFC65D","white",'gray33'))(6))
p_deduper <- plot_grid(p_D6, p_BC, ncol = 1, align = "v")
p_deduper

topptx(p_deduper,"GM_deduper.pptx",width = 11,height = 7)

dedup <- deduper_file %>%
  filter(method == "GM")
min(dedup$D6_rmse, na.rm = TRUE)
max(dedup$D6_rmse, na.rm = TRUE)
min(dedup$BC_rmse, na.rm = TRUE)
max(dedup$BC_rmse, na.rm = TRUE)

# calibrator
plot_GM_calibrator <- function(calibrator_file, metrics, breaks, color){
  calibrator <- calibrator_file %>%
    filter(method == "GM") %>%
    mutate(tools_aligner_counter_dedup_cali = paste(tools_aligner, tools_counter, tools_deduper, tools_calibrator, sep = "_")) %>%
    select(tools_aligner_counter_dedup_cali, metrics) %>%
    as.data.frame()
  rownames(calibrator) <- calibrator$tools_aligner_counter_dedup_cali
  calibrator <- calibrator[-1]
  
  p <- pheatmap(calibrator,
                cluster_cols = F,
                cluster_rows = F,
                annotation_legend = T,
                treeheight_row = 6,treeheight_col = 7,
                scale = "none",
                # scale = "column",
                # gaps_row = c(12),
                # gaps_col = c(24, 42),
                breaks = breaks,
                color = color,
                cellwidth = 14, cellheight = 14,
                clustering_distance_cols="correlation",
                clustering_distance_rows="correlation",
                #clustering_distance_cols = "euclidean",
                #clustering_distance_rows = "euclidean",
                show_colnames =F,show_rownames = T,
                #annotation_row=anno_20,
                #annotation_col=anno_counter,
                #  border=F,
                border_color = "black")
  g <- grid::grid.grabExpr(print(p))
  return(g)
}

p_D6 <- plot_GM_calibrator(calibrator_file, "D6_rmse", breaks = seq(13.8, 14.2, by = 0.1), color = colorRampPalette(c("#FFC65D","white",'gray33'))(5))
p_BC <- plot_GM_calibrator(calibrator_file, "BC_rmse", breaks = seq(14.1, 14.6, by = 0.1), color = colorRampPalette(c("#FFC65D","white",'gray33'))(6))
p_calibrator <- plot_grid(p_D6, p_BC, ncol = 1, align = "v")
p_calibrator

topptx(p_calibrator,"GM_calibrator.pptx",width = 11,height = 7)

cali <-  calibrator_file %>%
  filter(method == "GM")
min(cali$D6_rmse, na.rm = TRUE)
max(cali$D6_rmse, na.rm = TRUE)
min(cali$BC_rmse, na.rm = TRUE)
max(cali$BC_rmse, na.rm = TRUE)
