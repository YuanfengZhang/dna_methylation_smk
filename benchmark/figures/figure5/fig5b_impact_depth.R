# 2025-0826
# EQA impact_depth
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(dplyr)
library(tidyr)
library(pheatmap)

depth <- read.table("beta_depth_rows.csv",sep = ",",header = T)
head(depth)
colnames(depth)[5] <- "Depth_bins"
colnames(depth)[6] <- "Beta_bins"
# three bins
depth <- depth %>%
  mutate(Depth_bins = case_when(
    Depth_bins %in% c("1x", "2x", "3x", "4x") ~ "1–5x",
    Depth_bins == "5–10x" ~ "5–10x",
    Depth_bins == "≥10x" ~ "≥10x"
  ))

depth$Depth_bins <-  factor(depth$Depth_bins, levels = c("1–5x", "5–10x","≥10x"))
depth$Beta_bins <-  factor(depth$Beta_bins, levels = c("0–10","10–20","20–30","30–40","40–50",
                                                     "50–60","60–70","70–80","80–90","90–100"))
dodge_width <- 0.9
p <- ggplot(depth, aes(x = Beta_bins, y = rmse, fill = Depth_bins)) +
  geom_violin(width = 1.5, alpha = 0.5, color = NA, position = position_dodge(width = dodge_width)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", alpha = 0.7,
               position = position_dodge(width = dodge_width)) +
  scale_fill_manual(values = c("#f7b293","#FFDDAD","#AEC5EB","#3A405A","#BEE5A0","#049A8F")) +
  coord_cartesian(ylim = c(0, 65)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  facet_grid2(~ Beta_bins, scales = "free_x", space = "free_x") +
  labs(x = "Beta Bins", y = "RMSE");p

topptx(p,"impact_depth.pptx",width = 11,height =3)
max(depth$rmse)

# six bins
depth$Depth_bins <-  factor(depth$Depth_bins, levels = c("1x", "2x", "3x", "4x", "5–10x","≥10x"))
depth$Beta_bins <-  factor(depth$Beta_bins, levels = c("0–10","10–20","20–30","30–40","40–50",
                                                       "50–60","60–70","70–80","80–90","90–100"))
dodge_width <- 0.9
p <- ggplot(depth, aes(x = Beta_bins, y = rmse, fill = Depth_bins)) +
  geom_violin(width = 1.5, alpha = 0.5, color = NA, position = position_dodge(width = dodge_width)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "black", alpha = 0.7,
               position = position_dodge(width = dodge_width)) +
  scale_fill_manual(values = c("#f7b293","#FFDDAD","#AEC5EB","#3A405A","#BEE5A0","#049A8F")) +
  coord_cartesian(ylim = c(0, 65)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  facet_grid2(~ Beta_bins, scales = "free_x", space = "free_x") +
  labs(x = "Beta Bins", y = "RMSE");p

topptx(p,"impact_depth2.pptx",width = 11,height =3)

# 2025-0914
# EQA impact_factors beta_depth
com <- read.table("beta_depth_rows_compact.csv",sep = ",",header = T)
head(com)

fig5b <-com %>%
  filter(feature == "compact_beta_bin:grouped_depth_bin") %>%
  separate(fgroup, into = c("beta", "depth"), sep = ":")

fig5b$depth <-  factor(fig5b$depth, levels = c("1–4x", "5–10x", "≥10x"))
fig5b$beta <- factor(fig5b$beta, levels = c("0–10","10–20","20–30","30–40","40–50",
                                            "50–60","60–70","70–80","80–90","90–100"))

fig5b <- fig5b %>%
  group_by(depth, beta) %>%
  summarise(mean = mean(rmse, na.rm=T),
            sd = sd(rmse, na.rm =T),
            .groups = "drop")

p <- ggplot(fig5b, aes(x= beta, y =mean, fill = depth)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c("#f7b293","#BEE5A0","#049A8F")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  facet_grid2(~ beta, scales = "free_x", space = "free_x") +
  labs(x = "Beta Bins", y = "RMSE");p

topptx(p,"impact_depth.pptx",width = 11,height =3)
