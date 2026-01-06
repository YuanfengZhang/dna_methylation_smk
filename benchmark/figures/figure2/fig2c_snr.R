# fig 3a
# SNR11
# 2024-1203

library(ggplot2)
library(ggh4x)
library(RColorBrewer)
library(eoffice)
snr <-read.csv("unified_snr.csv",sep = ",")
snr11 <- snr %>%
  mutate(number = log10(length + 1)) %>%
  pivot_longer(
  cols = c(SNR, number),
  names_to = "Metrics",
  values_to = "Values"
)
snr11$Platform <- sapply(strsplit(as.character(snr11$lab),""),function(x){paste(x[1], x[2], sep ="")})
snr11$Platform <- gsub("BS","WGBS",snr11$Platform)
snr11$Platform <- gsub("EM","EM-seq",snr11$Platform)
snr11$Platform <- gsub("GM","GM-seq",snr11$Platform)
snr11$Platform <- gsub("RR", "RRBS",snr11$Platform)
snr11$Platform <- gsub("RM", "RRMS",snr11$Platform)
snr11$Platform <- factor(levels = c("WGBS", "EM-seq", "RRBS", "GM-seq", "RRMS", "MA"),snr11$Platform)
snr11$Metrics <- factor(snr11$Metrics,levels = c("SNR","number"))
snr11$lab <- factor(levels =c("BS1", "BS2", "BS3", "BS4",
                              "EM1", "EM2", "EM3", "EM4",
                              "RR1",
                              "GM1", "GM2", "GM3",
                              "RM1",
                              "MA1", "MA2", "MA3"), snr11$lab)
line_values <- snr11 %>%
  group_by(Metrics) %>%
  summarise(line_value = mean(Values)-sd(Values), .groups = "drop")  
method_colors <- c("WGBS" = "#FFDDAD",
                   "EM-seq" = "#AEC5EB",
                   "RRBS" = "#049A8F",
                   "GM-seq" = "#3A405A",
                   "RRMS" = "#BEE5A0",
                   "MA"= "#3C5487")

p <- ggplot(snr11, aes(x=lab, y=Values, fill=lab)) +
  geom_bar(stat = "identity",color="black", position = position_dodge(width = 0.8),
           linewidth=0.5 ,width = 0.8) +
  geom_text(aes(label = sprintf("%.1f", Values), y = Values), position = position_dodge(width = 0.3),
            size = 3,
            vjust = -1) +  
  facet_grid2(Metrics~Platform,scales = "free",space="free_x",independent = "y",
              strip  = strip_nested(
                background_x = elem_list_rect(fill = method_colors,color = NA),
                text_x = element_text(color ="black",face = "bold"),
                background_y = elem_list_rect(fill =  
                                                c('white'),color = "white"),
                text_y = element_text(color ="black",face = "bold")))+
  geom_hline(data=line_values,
             aes(yintercept=line_value),
             linewidth=0.5, colour="#C0504D", lty=5)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold",size = 18),
        axis.text.x = element_text(face = "bold", color = "black",size = 14),
        axis.text.y = element_text(color = "black", face = "bold",size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", color = "black",size = 14),
        axis.ticks.x =  element_blank(),
        legend.title = element_text(color = "black", face = "bold"),
        legend.text = element_text(color = "black", face = "bold"),
        legend.position = "none",
        legend.background = element_blank(),
        #panel.grid = element_line(colour = NA),
        panel.grid.minor  =  element_blank(),
        panel.grid.major.x  =  element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        #   panel.background = element_rect(fill = c("gray","blue"), colour = NA),
        strip.background.x = element_rect(fill = "white", colour = "white"),
        strip.text.x = element_text(colour = "black", face = "bold",size = 12),
        strip.text.y = element_text(colour = "black", face = "bold",size = 12),
        strip.placement = "inside",
        strip.switch.pad.grid = unit(1, "inch")) +
  scale_fill_viridis_d(option = "C", direction = -1)+
  labs(x = "Lab", y = "SNR")+
  facetted_pos_scales(
    y = list(
      Metrics == "SNR" ~ scale_y_continuous(limits = c(0, 38)),
      Metrics == "number" ~ scale_y_continuous(limits = c(0, 9))))+
  force_panelsizes(rows = c(2, 1));p

topptx(p,"snr_number.pptx", width = 11, height =7)
ggsave(p, filename = "Fig2c_SNR11_1x_bar.pdf",dpi=300, width =7, height = 3.5)