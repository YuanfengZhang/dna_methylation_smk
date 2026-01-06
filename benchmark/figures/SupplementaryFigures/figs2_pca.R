# SNR of each batch, depth=1x
library(dplyr)
library(ggplot2)
library(patchwork)

pcs <- read.csv("unified_pcs.csv")
snr <- read.csv("unified_snr.csv")
print(unique(pcs$lab))
pcs$lab <- factor(pcs$lab, levels = c("BS1","BS2","BS3","BS4",
                                      "EM1","EM2","EM3","EM4",
                                      "RR1","GM1","GM2","GM3",
                                      "RM1","MA1","MA2","MA3"))

type_colors <- c("D5"='#4CC3D9',"D6"='#7BC8A4',"F7"='#FFC65D',"M8"='#F16745',
                 "T1"="gray100","T2"="gray66","T3"="gray33","T4"="gray0")

plot_pca <- function(pcs, snr) {
  labs <- levels(pcs$lab)[levels(pcs$lab) %in% unique(pcs$lab)]
  plots <- list()
  # global_min_rmse <- min(data$rmse, na.rm = TRUE)
  # global_max_rmse <- max(data$rmse, na.rm = TRUE)
  
  for (i in seq_along(labs)) {
    lab <- labs[i]
    pcs_lab <- pcs %>% filter(lab == {{lab}}) %>%
      mutate(type = substr(sample,1,2))
    
    snr_lab <- snr %>% filter(lab == {{lab}})
    SNR_n <- round(snr_lab$SNR, 1)
    lens <- format(snr_lab$length, big.mark = ",", scientific = FALSE)

    p <- ggplot(pcs_lab, aes(x=PC1, y=PC2,
                             fill=type)) +
      geom_point(shape=21,size=2) +
      
      theme_bw() +
      theme(
        axis.text = element_text(),
        legend.position = "right",
        legend.background = element_blank(),
        plot.title = element_text(hjust=0.5, color="black", face="bold"),
        plot.subtitle = element_text(hjust=0.5)
      ) +
      guides(
        fill = guide_legend("Type", override.aes = list(shape=21, fill=type_colors))
      ) +
      scale_fill_manual("Type", values=type_colors) +
      labs(
        x=sprintf("PC1 (%.2f%%)", snr_lab$PC1*100),
        y=sprintf("PC2 (%.2f%%)", snr_lab$PC2*100),
        title = paste0(lab),
        subtitle=paste0("SNR=", SNR_n," (",lens,")")
      ) +
      theme(aspect.ratio=1/1)
    
    plots[[i]] <- p  
  }
  final_plot <- wrap_plots(plots, ncol = 4) +
    plot_layout(guides = "collect")
  
  return(final_plot)
}

p <- plot_pca(pcs, snr)
print(p)  

topptx(p,"FigS2.pptx",width = 12,height = 12) 
