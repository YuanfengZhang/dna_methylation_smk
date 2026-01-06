# 2025-0827
# EQA methcali_rmse
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(dplyr)
library(tidyr)
library(pheatmap)

before_rmse <- read.table("before_rmse.csv",sep = ",",header = T)
after_rmse <- read.table("after_rmse.csv",sep = ",",header = T)
head(rmse_beta)
rmse <- merge(before_rmse, after_rmse, 
              by = c("lab", "sample", "feature", "fgroup"),
              all = TRUE) %>%
  as_tibble() %>% 
  rename(rmse_original = rmse.x,
         rmse_cali = rmse.y,
         count = count.y) %>%
  select(-count.x)

# plot rmse_beta_bins
rmse_beta <- rmse %>%
  filter(feature == "compact_beta_bin") %>%
  pivot_longer(cols = c("rmse_original","rmse_cali"),
               names_to = "method",
               values_to ="Values")

rmse_beta$method <-  factor(rmse_beta$method, levels = c("rmse_original","rmse_cali"))
rmse_beta$fgroup <-  factor(rmse_beta$fgroup, levels = c("0–10","10–20","20–30","30–40","40–50",
                                                     "50–60","60–70","70–80","80–90","90–100"))
dodge_width <- 1
p <- ggplot(rmse_beta, aes(x = fgroup, y = Values, fill = method)) +
  geom_violin(width = 1, alpha = 0.2, color = NA, position = position_dodge(width = dodge_width)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", alpha = 0.5,
               position = position_dodge(width = dodge_width)) +
  scale_fill_manual(values = c("#0b527d","#de4649")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  facet_grid2(~ fgroup, scales = "free_x", space = "free_x") +
  labs(x = "Beta Bins", y = "RMSE");p

topptx(p,"methcali_rmse_beta.pptx",width = 8,height =3)

# plot rmse_depth_bins
rmse_depth <- rmse %>%
  filter(feature == "compact_depth_bin") %>%
  pivot_longer(cols = c("rmse_original","rmse_cali"),
               names_to = "method",
               values_to ="Values")

rmse_depth$method <-  factor(rmse_depth$method, levels = c("rmse_original","rmse_cali"))
rmse_depth$fgroup <-  factor(rmse_depth$fgroup, levels = c("1x", "2x", "3x", "4x", "5–10x","≥10x"))
dodge_width <- 1
p <- ggplot(rmse_depth, aes(x = fgroup, y = Values, fill = method)) +
  geom_violin(width = 1, alpha = 0.2, color = NA, position = position_dodge(width = dodge_width)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", alpha = 0.5,
               position = position_dodge(width = dodge_width)) +
  scale_fill_manual(values = c("#0b527d","#de4649")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 1),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "grey90", size = 0.3)
  ) +
  facet_grid2(~ fgroup, scales = "free_x", space = "free_x") +
  labs(x = "Depth Bins", y = "RMSE");p

topptx(p,"methcali_rmse_depth.pptx",width = 8,height =3)

# 0919 methcali rmse_depth
head(dat1)
dat1 <- RMSE_before[RMSE_before$feature=="compact_depth_bin",]
dat1$Type <- factor(dat1$Type ,levels=c("Before","After"))
dat1$fgroup <- factor(dat1$fgroup, levels = c("1x", "2x", "3x", "4x", "5–10x","≥10x"))

p<-ggplot(data=dat1,aes(x=fgroup,
                        y=rmse))+
  geom_boxplot(aes(fill=Type),
               outlier.shape = NA,
               outlier.colour = "gray75",
               outlier.alpha = 0.2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text =element_text(color = "black"),
        axis.text.x = element_text(angle=0,hjust=0.5,vjust=1),
        legend.position = "bottom")+
  scale_fill_manual(values = c("Before"="#2f81b7",
                               "After"="#c9211a"),
                    # c("Before"="#7670AA",
                    #  "After"="#D56"),
                    name="Group")+
  labs(x=NULL,y="RMSE (%)")+
  scale_y_continuous() +
  facet_grid(~fgroup, scales="free_x");p

# 计算 P 值
pvalue.df <- dat1 %>%
  group_by(fgroup) %>%
  summarise(pvalue = wilcox.test(rmse ~ Type)$p.value,
            .groups = "drop") %>%
  mutate(new_p = ifelse(pvalue == 0, min(pvalue[pvalue!=0]), pvalue))

# P 值图
p.top <- ggplot(pvalue.df, aes(x=fgroup, y=1, color=-log10(new_p))) +
  geom_point(size=5, shape=15) +
  scale_color_distiller(palette="Greys", direction=1, name="-log10(Pval)") +
  theme_void() +
  facet_grid(~fgroup, scales="free_x", space="free_x") +
  theme(legend.position="top",
        legend.key.size = unit(5,'mm'),
        legend.text = element_text(size=10)) +
  guides(color=guide_colorbar(title.position="top",
                              title.hjust=0.5,
                              barwidth=20));p.top

# 上下组合
p_all <- p.top + p+ plot_layout(ncol=1, heights=c(1,5))
p_all

topptx(p_all,"methcali_rmse_depth.pptx",height = 6,width = 8)


# 0918 methcali rmse_sample
library(ggpubr)
RMSE_before <- read.csv("before_rmse.csv")
RMSE_before$Type <-"Before"
After_before <- read.csv("after_rmse.csv")
After_before$Type <-"After"
RMSE_before <- rbind(RMSE_before,After_before)
head(RMSE_before)
table(RMSE_before$feature)
#RMSE_before$feature=="compact_beta_bin"
dat1 <- RMSE_before[RMSE_before$fgroup=="≥10x",]
dat1 <- dat1 %>%
  mutate(sample = substr(sample,1,2))
head(dat1)
dat1$Type <- factor(dat1$Type ,levels=c("Before","After"))
dat1 <- dat1 %>%
  mutate(Group = case_when(
    sample %in% c("D5","D6","F7","M8") ~ "Quartet sample",
    sample %in% c("T1","T2","T3","T4") ~ "tiltrated mixture samples",
    sample %in% c("BC","BL")           ~ "tumor-normal"
  ))
dat1$sample <- factor(dat1$sample, levels = unique(dat1$sample))
p<-ggplot(data=dat1,aes(x=sample,
                        y=rmse))+
  geom_boxplot(aes(fill=Type),
               outlier.shape = NA,
               outlier.colour = "gray75",
               outlier.alpha = 0.2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text =element_text(color = "black"),
        axis.text.x = element_text(angle=0,hjust=0.5,vjust=1),
        legend.position = "bottom")+
  scale_fill_manual(values = c("Before"="#2f81b7",
                               "After"="#c9211a"),
                    # c("Before"="#7670AA",
                    #  "After"="#D56"),
                    name="Group")+
  labs(x=NULL,y="RMSE (%)")+
  scale_y_continuous() +
  facet_grid(~Group, scales="free_x");p

# 计算 P 值
pvalue.df <- dat1 %>%
  group_by(Group, sample) %>%
  summarise(pvalue = wilcox.test(rmse ~ Type)$p.value,
            .groups = "drop") %>%
  mutate(new_p = ifelse(pvalue == 0, min(pvalue[pvalue!=0]), pvalue))

# P 值图
p.top <- ggplot(pvalue.df, aes(x=sample, y=1, color=-log10(new_p))) +
  geom_point(size=5, shape=15) +
  scale_color_distiller(palette="Greys", direction=1, name="-log10(Pval)") +
  theme_void() +
  facet_grid(~Group, scales="free_x", space="free_x") +
  theme(legend.position="top",
        legend.key.size = unit(5,'mm'),
        legend.text = element_text(size=10)) +
  guides(color=guide_colorbar(title.position="top",
                              title.hjust=0.5,
                              barwidth=20))

# 上下组合
p_all <- p.top + p+ plot_layout(ncol=1, heights=c(1,5))
p_all

topptx(p_all,"methcali_rmse_sample.pptx",height = 6,width = 8)
