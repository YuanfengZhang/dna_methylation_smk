# 2025-0914
# fig6b
library(ggpubr)
RMSE_before <- read.csv("before_rmse.csv")
RMSE_before$Type <-"Before"
After_before <- read.csv("after_rmse.csv")
After_before$Type <-"After"
RMSE_before <- rbind(RMSE_before,After_before)
head(RMSE_before)
table(RMSE_before$feature)
#RMSE_before$feature=="compact_beta_bin"
dat1 <- RMSE_before[RMSE_before$feature=="compact_beta_bin",]
head(dat1)
dat1$Type <- factor(dat1$Type ,levels=c("Before","After"))

p<-ggplot(data=dat1,aes(x=fgroup,
                        y=rmse))+
  geom_boxplot(aes(fill=Type),
               outlier.shape = 1,
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
  scale_y_continuous(breaks = c(0,10,20,30),
                     labels = c("0%","10%","20%","30%"));p
# error bar---
ggplot_build(p)$data[[1]] %>% select(x,ymin,ymax) -> errorbar.df

p.bottom<-p+
  geom_segment(data = errorbar.df,
               aes(x=x-0.15,xend=x+0.15,y=ymin,yend=ymin))+
  geom_segment(data = errorbar.df,
               aes(x=x-0.15,xend=x+0.15,y=ymax,yend=ymax));p.bottom
# error bar---

dat1 %>% 
  pull(fgroup) %>% 
  unique() -> group.info

pvalue.df<-tibble(x=character(),
                  pvalue=numeric())
for(info in group.info){
  dat1 %>% 
    filter(fgroup==info) -> tmp.df
  wilcox.test(rmse~Type,
              data=tmp.df) -> a
  add_row(pvalue.df,
          x=info,
          pvalue=a$p.value) -> pvalue.df
}

min_p<-min(pvalue.df %>% filter(pvalue !=0) %>% 
             pull(pvalue))
pvalue.df %>% 
  mutate(new_p=case_when(
    pvalue == 0 ~ min_p,
    TRUE ~ pvalue
  )) ->pvalue.df

ggplot_build(p)$data[[1]] %>% select(x,ymin,ymax) -> errorbar.df

# pvalue

library(RColorBrewer)
scale_color_distiller()
p.top<-ggplot(data = pvalue.df,aes(x=x,y=1))+
  geom_point(size=5,shape=15,
             aes(color=-log10(new_p)))+
  scale_color_distiller(palette = "Greys",
                        direction = 1,
                        name="-log10(Pval)")+
  theme_void()+
  theme(legend.position = "top",
        legend.key.size = unit(5,'mm'),
        legend.text = element_text(size=10))+
  guides(color=guide_colorbar(title.position = "top",
                              title.hjust = 0.5,
                              barwidth = 20));p.top
# combine
library(patchwork)

p_all <- p.top+
  p.bottom+
  plot_layout(ncol = 1,heights = c(1,5));p_all
p_all <- as_ggplot(p_all)
topptx(p_all,"Fig6_b_0914_3.pptx",height = 5,width = 10)
ggsave(p_all,"Fig6_b_0914_3.png",height = 7,width = 12)
