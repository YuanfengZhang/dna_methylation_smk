# 2025-0825
# EQA impact_factors
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(dplyr)
library(tidyr)
library(pheatmap)

fac <- read.table("impact_factors.csv",sep = ",",header = T)
head(fac)

fac$feature[fac$feature == ""] <- "BWT Compression Ratio"
fac$feature_group <- gsub("sequencing depth", "Depth", fac$feature_group)
fac$feature_group <- gsub("CpG density", "CpG", fac$feature_group)
fac$feature_group <- gsub("detailed GC content", "GC content", fac$feature_group)
fac$feature_group <- gsub("sequence complexity", "Sequence", fac$feature_group)

fac$feature <- factor(fac$feature, levels = c("Sequencing Depth",
                                              "CpG Density",
                                              "b1","b2","b3","b4","b5",
                                              "a2","a3","a4","a5",
                                              "GC%","GC skew","CpG / GC",
                                              "Shannon Entropy","BWT Compression Ratio",
                                              "Enhancer","Promoter","Genetic Location"
                                              ))
fac$feature_group <- factor(fac$feature_group, levels = c("Depth",
                                                          "CpG",
                                                          "motif",
                                                          "GC content",
                                                          "Sequence",
                                                          "genomic locations"
                                                          ))
fac <- fac %>%
  mutate(importance = log10(importance + 0.01)) 
"#ecdf42"
"#efc47a"
"#e8a3c1"
"#968abb"
"#b8a6ca"
"#a1c7cd"
"#7dc6e3"
"#80b778"

p <- ggplot(data = fac, aes(x = feature, y = importance, fill = feature_group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("#efc47a",
                               "#e8a3c1",
                               "#b8a6ca",
                               "#a1c7cd",
                               "#7dc6e3",
                               "#80b778")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    axis.title  = element_text(color = "black"),
    strip.text = element_text(),
    panel.grid.minor  = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) +
  facet_grid2(~feature_group, scales = "free",space="free_x",
              strip  = strip_nested(
                background_x = elem_list_rect(fill =  
                                                c("#efc47a",
                                                  "#e8a3c1",
                                                  "#b8a6ca",
                                                  "#a1c7cd",
                                                  "#7dc6e3",
                                                  "#80b778"),color = NA))) +
  
  labs(y = "Importance", x = "Feature")+ 
  scale_y_continuous(limits = c(-2.3, 0.8));p

topptx(p,"factors.pptx",width = 10,height =5)

# RR1
RR1 <- fac %>%
  filter(lab == "RR1")
head(RR1)

p <- ggplot(data = RR1, aes(x = feature, y = importance, fill = feature_group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("#f7b293","#FFDDAD","#AEC5EB","#3A405A","#BEE5A0","#049A8F")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    axis.title  = element_text(color = "black"),
    strip.text = element_text(),
    panel.grid.minor  = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  ) +
  facet_grid2(~feature_group, scales = "free",space="free_x",independent = "y",
              strip  = strip_nested(
                background_x = elem_list_rect(fill =  
                                                c("#f7b293","#FFDDAD","#AEC5EB","#3A405A","#BEE5A0","#049A8F"),color = NA))) +
  labs(y = "Importance", x = "Feature");p

topptx(p,"RR1_factors.pptx",width = 10,height =5)
