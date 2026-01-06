library(tidyverse)
library(glue)
library(eoffice)
unified_dir <- "."
spearmanr <- read.csv(glue("{unified_dir}/unified_spearmanr.csv"))

# 1. Define the lab order
lab_order <- c(
    "BS1", "BS2", "BS3", "BS4",
    "EM1", "EM2", "EM3", "EM4",
    "RR1", "GM1", "GM2", "GM3",
    "MA1", "MA2", "MA3", "RM1"
)
# 1. Define the label order
label_order <- c(
    "BL", "BC",
    "T4", "T3", "T2", "T1",
    "M8", "F7", "D6", "D5"
)

spearmanr <- spearmanr %>% filter(label != "HF")
spearmanr$lab <- factor(spearmanr$lab, levels = lab_order)
spearmanr$label <- factor(spearmanr$label, levels = label_order)

min_spr <- min(spearmanr$spearmanr, na.rm = TRUE)
max_spr <- max(spearmanr$spearmanr, na.rm = TRUE)

spr <- ggplot(spearmanr, aes(x = lab, y = label)) +
    # Create circles with color based on spearmanr values
    geom_point(aes(fill = spearmanr), shape = 21, size = 10, stroke = NA) +
    # Add text labels with 2 decimal places
    geom_text(aes(label = sprintf("%.2f", spearmanr), 
                  color = ifelse(spearmanr >= 0.75,
                                 "white",
                                 "black")),
              size = 3.5) +
    scale_color_manual(values = c("white" = "white",
                                  "black" = "black")) +
    # Set color gradient from white (0) to green (1)
    scale_fill_gradient(low = "white", high = "#509322",
                        limits = c(min_spr, max_spr)) +
    # 2. Move x-axis labels to the top
    scale_x_discrete(position = "top") +
    # 3. Increase spacing between rows
    scale_y_discrete(expand = expansion(add = 0.7)) +
    # Customize appearance
    theme_minimal() +
    theme(
        panel.grid = element_blank()
    ) +
    labs(
        x = "Lab", 
        y = "label",
        fill = "SpearmanR"
    )

topptx(spr, glue("{unified_dir}/unified_spearmanr_circles.pptx"),
       width = 10, height = 4)

jaccard <- read.csv(glue("{unified_dir}/unified_jaccard.csv"))
jaccard <- jaccard %>% filter(label != "HF")
jaccard$lab <- factor(jaccard$lab, levels = lab_order)
jaccard$label <- factor(jaccard$label, levels = label_order)
min_jac <- min(jaccard$jaccard, na.rm = TRUE)
max_jac <- max(jaccard$jaccard, na.rm = TRUE)
jac <- ggplot(jaccard, aes(x = lab, y = label)) +
    # Create circles with color based on jaccard values
    geom_point(aes(fill = jaccard), shape = 21, size = 10, stroke = NA) +
    # Add text labels with 2 decimal places
    geom_text(aes(label = sprintf("%.2f", jaccard)),
              size = 3.5) +
    # Set color gradient from white (0) to green (1)
    scale_fill_gradient(low = "white", high = "#C0392B",
                        limits = c(min_jac, max_jac)) +
    # 2. Move x-axis labels to the top
    scale_x_discrete(position = "top") +
    # Customize appearance
    theme_minimal() +
    theme(
        panel.grid = element_blank()
    ) +
    labs(
        x = "Lab", 
        y = "label",
        fill = "Jaccard"
    )

topptx(jac, glue("{unified_dir}/unified_jaccard_circles.pptx"),
       width = 10, height = 4)