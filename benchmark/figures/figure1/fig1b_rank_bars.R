library(tidyverse)
library(glue)
library(eoffice)
library(extrafont)
loadfonts(device = "win")

setwd(".")
accuracy_dir <- "."


rmse <- data.frame(
    rank = 1:17,
    lab = c("MA2", "MA1", "MA3",
            "GM3", "RR1", "GM1",
            "BS1", "GM2", "BS3",
            "EM2", "EM4", "EM3",
            "RM1", "BS2", "EM1",
            "NP1", "BS4"),
    penalty = c(3, 17, 42, 63,  100,
                109, 114, 120, 168,
                198, 205, 208, 235,
                242, 276, 300, 320)
)

spearmanr <- data.frame(
    rank = 1:15,
    lab = c("MA2", "MA1",
            "RR1", "GM3", "GM1",
            "GM2", "EM3", "BS2",
            "EM2", "EM4", "BS3",
            "RM1", "EM1", "BS4", "BS1"),
    penalty = c(14, 19, 39, 47,
                63, 71, 75, 91, 98,
                114, 127, 131, 135,
                144, 152)
)

rmse$lab <- factor(rmse$lab,
                   levels = rev(rmse$lab[order(rmse$rank)]))
spearmanr$lab <- factor(spearmanr$lab,
                        levels = rev(spearmanr$lab[order(spearmanr$rank)]))

rmse_bars <- ggplot(rmse, aes(x = penalty, y = lab)) +
    geom_col(width = 0.8, fill = "grey") +
    geom_text(aes(label = penalty), 
              hjust = -0.1, 
              size = 3.5) +
    scale_x_log10(
        name = "Penalty (log scale)",
        breaks = c(1, 3, 10, 30, 100, 300),
        labels = c(1, 3, 10, 30, 100, 300)
    ) +
    labs(y = "Lab", title = "Penalty by Lab") +
    theme_minimal() +
    theme(
        text = element_text(family = "Arial"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()
    )

spearmanr_bars <- ggplot(spearmanr,
                         aes(x = penalty, y = lab)) +
    geom_col(width = 0.8, fill = "grey") +
    geom_text(aes(label = penalty),
              hjust = -0.1,
              size = 3.5) +
    scale_x_log10(
        name = "Penalty (log scale)",
        breaks = c(10, 20, 50, 100, 200),
        labels = c(10, 20, 50, 100, 200)
    ) +
    labs(y = "Lab", title = "Penalty by Lab (SpearmanR)") +
    theme_minimal() +
    theme(
        text = element_text(family = "Arial"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()
)

topptx(rmse_bars,
       glue("{accuracy_dir}/rmse_rank_bars.pptx"),
       width = 2, height = 4)
topptx(spearmanr_bars,
       glue("{accuracy_dir}/spearmanr_rank_bars.pptx"),
       width = 2, height = 4)
