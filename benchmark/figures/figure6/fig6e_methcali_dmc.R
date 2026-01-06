# 2025-0828
# EQA methcali_dmc
install.packages("ggdist")
if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github('jorvlan/raincloudplots')

library(raincloudplots)
library(ggrain)
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(dplyr)
library(tidyr)
library(pheatmap)

dmc <- read.table("dmc_merged.csv",sep = ",",header = T)
head(dmc)
dmc <-  dmc %>%
  filter(p_type == "p",
         tool %in% c('methylkit', 'methylsig.beta_binomial', 'methylsig.binomial', 'methylsig.dss'),
         ! is.na(mcc))
print(unique(dmc$tool))

dmc$sample_pair <- factor(dmc$sample_pair,levels = c("BC_vs_BL", "D6_vs_D5", "D6_vs_F7", "D6_vs_M8"))
dmc$tool <- factor(dmc$tool,levels = c("methylkit", "methylsig.beta_binomial", "methylsig.binomial", "methylsig.dss"))
dmc$treatment <- factor(dmc$treatment,levels = c("before", "after"))

dodge_width <- 1
p <- ggplot(data =dmc,       
            aes(x=sample_pair, y=mcc, fill=treatment)) +
  geom_violin(width = 0.8, alpha = 0.2, color = NA, position = position_dodge(width = dodge_width)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", alpha = 0.6,
               position = position_dodge(width = dodge_width)) +
  scale_fill_manual(values = c("#0b527d","#de4649"))+
  theme_bw()+
  theme(
    axis.text.x = element_text(face = "bold",color = "black"),
    axis.text.y = element_text(face = "bold",color = "black"),
    axis.title  = element_text(face = "bold",color = "black"),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    panel.grid.minor  =  element_blank(),
    panel.spacing = unit(0.2, "cm"),
    legend.background =  element_blank(),
    legend.title = element_text(face = "bold",color = "black"),
    legend.text = element_text(color = "black")
  )+
  facet_grid2(.~tool, scales = "free",space="free_x",independent = "y",
              strip  = strip_nested(
                background_x = elem_list_rect(
                  fill = "white",color = "white"),
                background_y = elem_list_rect(
                  fill = "white",color ="white")))+
  labs(y="MCC",x="Sample_pair")+
  facetted_pos_scales(
    y = list(
      tool == "methylkit" ~ scale_y_continuous(),
      tool == "methylsig.beta_binomial" ~ scale_y_continuous(),
      tool == "methylsig.binomial" ~ scale_y_continuous(),
      tool == "methylsig.dss" ~ scale_y_continuous()
      ));p

topptx(p,"dmc_cali.pptx",width = 8,height =3)

# plot rainclound
dmc_wider <- dmc %>%
  select("sample_pair", "lab","tool","treatment","mcc") %>%
  distinct(sample_pair, lab, tool, treatment, .keep_all = TRUE) %>%
  pivot_wider(names_from = treatment,
              values_from = mcc)
head(dmc_wider)

Mentored_pre_post <- data_1x1( 
  array_1 = dmc_wider$before, #first set of values
  array_2 = dmc_wider$after, #second set of values
  jit_distance = .09,
  jit_seed = 321) 

p <- raincloud_1x1_repmes(
  data = Mentored_pre_post,
  colors = (c("#0b527d","#de4649")), 
  fills = (c("#0b527d","#de4649")), 
  line_color = 'gray',
  line_alpha = .1,
  size = .5,
  alpha = 0.4,
  align_clouds = FALSE) +
  scale_x_continuous(breaks=c(1,2), labels=c("before", "after"), limits=c(0, 3)) +
  xlab("Method") + 
  ylab("RMSE") +
  theme_classic()
p

topptx(p,"dmc_cali2.pptx",width = 3,height =2)
