# 2025-0826
# EQA impact_den
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(dplyr)
library(tidyr)
library(pheatmap)

den <- read.table("cpg_density.csv",sep = ",",header = T)
head(den)
colnames(den)[7] <- "Density_bins"
colnames(den)[6] <- "Beta_bins"

den$Density_bins <-  factor(den$Density_bins, levels = c("CpG Islands", "CpG Shores","CpG Shelves","Open Sea"))
den$Beta_bins <-  factor(den$Beta_bins, levels = c("0–10","10–20","20–30","30–40","40–50",
                                                       "50–60","60–70","70–80","80–90","90–100"))
dodge_width <- 0.9
p <- ggplot(den, aes(x = Density_bins, y = rmse, fill = Density_bins)) +
  geom_violin(width = 1, alpha = 0.5, color = NA, position = position_dodge(width = dodge_width)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", alpha = 0.7,
               position = position_dodge(width = dodge_width)) +
  scale_fill_manual(values = c("#f7b293","#FFDDAD","#AEC5EB","#3A405A","#BEE5A0","#049A8F")) +
  #coord_cartesian(ylim = c(0, 65)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  # facet_grid2(~ Density_bins, scales = "free_x", space = "free_x") +
  labs(x = "Beta Bins", y = "RMSE");p

topptx(p,"impact_den.pptx",width = 6,height =3)
max(den$rmse)
