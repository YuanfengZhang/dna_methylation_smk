# 0821
# Cytosines in depth
library(dplyr)
library(ggplot2)
library(ggh4x)
library(eoffice)

cytosines <- read.csv("dce_stat.csv",sep = ",")
cytosines <- cytosines %>%
  rename("[0%, 5%)" = "X.5",
         "[5%, 10%)" = "X5_10",
         "≥10%" = "X..10") 
fig_c <- pivot_longer(
  cytosines,
  cols = c("[0%, 5%)", "[5%, 10%)","≥10%"),
  names_to = "Metrics",
  values_to = "Values"
)  
fig_c$Platform <- sapply(strsplit(as.character(fig_c$lab),""),function(x){paste(x[1], x[2], sep ="")})
fig_c$Platform <- gsub("BS","WGBS",fig_c$Platform)
fig_c$Platform <- gsub("EM","EM-seq",fig_c$Platform)
fig_c$Platform <- gsub("GM","GM-seq",fig_c$Platform)
fig_c$Platform <- gsub("RR", "RRBS",fig_c$Platform)
fig_c$Platform <- gsub("RM", "RRMS",fig_c$Platform)
fig_c$Platform <- factor(levels = c("WGBS", "EM-seq", "RRBS", "GM-seq", "RRMS", "MA"),fig_c$Platform)
fig_c$Metrics <- factor(fig_c$Metrics,levels = c("≥10%", "[5%, 10%)","[0%, 5%)"))
fig_c$lab <- factor(levels =c("BS1", "BS2", "BS3", "BS4",
                              "EM1", "EM2", "EM3", "EM4",
                              "RR1",
                              "GM1", "GM2", "GM3",
                              "RM1",
                              "MA1", "MA2", "MA3"), fig_c$lab)

cytosines$variable <-factor(levels = c("≥10%", "[5%, 10%)","[0%, 5%)"),cytosines$variable)

p <- ggplot(fig_c,aes(x=lab,y=Values,fill=Metrics))+  
  geom_bar(stat ="identity",alpha=0.9)+  
  facet_grid2(.~Platform,scales = "free",space="free_x",independent = "y",
              strip  = strip_nested(
                background_x = elem_list_rect(fill =  
                                                c('#f7dd82' ,"#f38181","#C5008499","#925E9FFF","#9900CCFF"),color = NA)))+
  theme_bw()+
  theme(axis.text.x =  element_text(color = "black",face = "bold"),
        axis.text.y = element_text(color = "black",face = "bold"),
        axis.title.x =  element_text(color = "black",face = "bold"),
        axis.title.y = element_text(color = "black",face = "bold"),
        axis.ticks.x = element_blank(),
        legend.title = element_text(color = "black",face = "bold"),
        legend.text = element_text(color = "black",face = "bold"),
        legend.position = "bottom",legend.background = element_blank(),
        panel.grid = element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA))+
  theme(strip.background.y = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(size = 12, colour = "black",face = "bold"),
        strip.text.y = element_text(size = 12, colour = "black",face = "bold")) + 
  theme(strip.placement = "inside") + 
  theme(strip.switch.pad.grid = unit(0.2, "inch"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  scale_fill_manual(values = c("#3C5488FF", "#8491B4FF", "#F39B7FFF"))+
  labs(
    fill="Depth distribution",
    y="Percentage (%)");p

topptx(p,"cytosines_depth.pptx", width = 11, height = 7)
ggsave(p, filename = "Fig2_b_MAD_bar.pdf",dpi=300, width =7, height = 3.5)
