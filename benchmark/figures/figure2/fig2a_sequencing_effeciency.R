# 2025-0818
# EQA fig2a-temp
# a--depth or cytosine efficiency
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
setwd(".")
de <- read.table("depth_efficiency.csv",sep = ",",header = T)
cy <- read.table("cytosine_efficiency.csv",sep = ",",header = T)
head(de)
head(cy)
de$Metrics <- "Depth"
cy$Metrics <- "Cytosine"
colnames(de)[7]<-"efficiency"
colnames(cy)[7]<-"efficiency"

fig2a <- rbind(de[,c(1:4,7:8)],cy[,c(1:4,7:8)])
head(fig2a)
table(fig2a$label)
fig2a$protocol <- sapply(strsplit(as.character(fig2a$lab),""),function(x){paste(x[1],x[2],sep = "")})
fig2a$label <- factor(fig2a$label,levels = c("D5","D6","F7","M8","BC","BL","T1","T2","T3","T4","HF"))
fig2a$Metrics <- factor(fig2a$Metrics,levels = c("Depth","Cytosine"))
fig2a$protocol <- gsub("BS","WGBS",fig2a$protocol)
fig2a$protocol <- gsub("EM","EMseq",fig2a$protocol)
fig2a$protocol <- gsub("GM","GM-seq",fig2a$protocol)
fig2a$protocol <- gsub("RR","RRBS",fig2a$protocol)
fig2a$protocol <- gsub("RM","RRMS",fig2a$protocol)

fig2a$protocol <- factor(fig2a$protocol,
                         levels = c("WGBS","EMseq","GM-seq","RRMS","RRBS"))

p <- ggplot(data =fig2a,       
            aes(x=lab, y=efficiency, fill=lab)) +
  geom_half_violin(side = "r", color=NA, alpha=0.35) +    
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, 
                    outlier.size = 1,
                    outlier.stroke = 0.2,
                    width=0.2, linewidth=0.5) +    
  geom_half_point_panel(aes( fill = label),side = "l", 
                        range_scale = .85,
                        shape=21, size=1.5, color="black")+
  scale_fill_manual(values = c('#4CC3D9' ,'#7BC8A4' ,'#FFC65D', '#F16745',"#69001f","#f7b293",
                               "gray100","gray66","gray33","gray0","#8039f9",
                               "#FFDDAD","#FFDDAD","#FFDDAD","#FFDDAD",
                               "#AEC5EB","#AEC5EB","#AEC5EB","#AEC5EB",
                               "#3A405A","#3A405A","#3A405A",
                               "#BEE5A0","#049A8F"))+
  theme_bw()+
  theme(
    axis.text.x = element_text(face = "bold",color = "black"),
    axis.text.y = element_text(face = "bold",color = "black"),
    axis.title  = element_text(face = "bold",color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor  =  element_blank(),
    legend.background =  element_blank(),
    legend.title = element_text(face = "bold",color = "black"),
    legend.text = element_text(color = "black")
  )+
  facet_grid2(Metrics~protocol,scales = "free",space="free_x",independent = "y",
              strip  = strip_nested(
                background_x = elem_list_rect(
                  fill = "white",color = "white"),
                background_y = elem_list_rect(
                  fill = "white",color ="white")))+
  labs(y="Efficiency",x="Assay")+
  facetted_pos_scales(
    y = list(
      Metrics == "Depth" ~ scale_y_continuous(labels = function(x) paste0(x, "X")),
      Metrics == "Cytosine" ~ scale_y_continuous(labels = scales::percent_format(scale = 1)) ));p


ggsave("Fig2a_cytosine_depth_eff.png",p,width = 12,height = 5.5)
topptx(p,"Fig2a_cytosine_depth_eff.pptx",width = 12,height = 5.5)

