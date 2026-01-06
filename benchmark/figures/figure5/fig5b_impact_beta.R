# 2025-0826
# EQA impact_beta
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(dplyr)
library(tidyr)
library(pheatmap)

beta <- read.table("beta_rows.csv",sep = ",",header = T)
head(beta)
colnames(beta)[3] <- "Beta_bins"
beta$method <- sapply(strsplit(as.character(beta$lab),""),function(x){paste(x[1],x[2],sep = "")})
beta$method <-  factor(beta$method, levels = c('BS', 'EM', 'RR', 'GM'))
beta$Beta_bins <- factor(beta$Beta_bins, levels = c("0–10","10–20","20–30","30–40","40–50",
                                                     "50–60","60–70","70–80","80–90","90–100"))
dodge_width <- 0.9
p <- ggplot(beta, aes(x = Beta_bins, y = RMSE, fill = method)) +
  geom_violin(width = 1.5, alpha = 0.5, color = NA, position = position_dodge(width = dodge_width)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", alpha = 0.7,
               position = position_dodge(width = dodge_width)) +
  scale_fill_manual(values = c("#f7b293","#FFDDAD","#AEC5EB","#3A405A","#BEE5A0","#049A8F")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  facet_grid2(~ Beta_bins, scales = "free_x", space = "free_x") +
  labs(x = "Beta Bins", y = "RMSE");p

topptx(p,"impact_beta.pptx",width = 11,height =3)

# 2025-0914
# EQA impact_factors beta
com <- read.table("beta_depth_rows_compact.csv",sep = ",",header = T)
head(com)

fig5b <-com %>%
  filter(feature == "compact_beta_bin:grouped_depth_bin") 
fig5b$method <- sapply(strsplit(as.character(fig5b$lab),""),function(x){paste(x[1],x[2],sep = "")})
fig5b$method <-  factor(fig5b$method, levels = c('BS', 'EM', 'RR', 'GM'))
fig5b$Beta_bins <- factor(fig5b$fgroup, levels = c("a","b","c","d","e"))

fig5b <- fig5b %>%
  group_by(method, Beta_bins) %>%
  summarise(mean = mean(rmse, na.rm=T),
            sd = sd(rmse, na.rm =T),
            .groups = "drop")

p <- ggplot(fig5b, aes(x= Beta_bins, y =mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c("#FFDDAD","#AEC5EB","#049A8F","#3A405A")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  facet_grid2(~ Beta_bins, scales = "free_x", space = "free_x") +
  labs(x = "Beta Bins", y = "RMSE");p

topptx(p,"impact_beta.pptx",width = 11,height =3)


# 0921 factor beta 10 bins
beta <- read.table("beta_rows.csv",sep = ",",header = T)
head(beta)
colnames(beta)[3] <- "Beta_bins"
beta$method <- sapply(strsplit(as.character(beta$lab),""),function(x){paste(x[1],x[2],sep = "")})
beta$method <-  factor(beta$method, levels = c('BS', 'EM', 'RR', 'GM'))
beta$Beta_bins <- factor(beta$Beta_bins, levels = c("0–10","10–20","20–30","30–40","40–50",
                                                    "50–60","60–70","70–80","80–90","90–100"))
head(beta)
beta <- beta %>%
  group_by(method, Beta_bins) %>%
  summarise(mean = mean(RMSE, na.rm=T),
            sd = sd(RMSE, na.rm =T),
            .groups = "drop")

p <- ggplot(beta, aes(x= Beta_bins, y =mean, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c("#FFDDAD","#AEC5EB","#049A8F","#3A405A")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  facet_grid2(~ Beta_bins, scales = "free_x", space = "free_x") +
  labs(x = "Beta Bins", y = "RMSE");p

topptx(p,"impact_beta.pptx",width = 11,height =3)
