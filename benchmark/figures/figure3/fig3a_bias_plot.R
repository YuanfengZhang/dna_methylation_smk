library(extrafont)
library(eoffice)
library(ggridges)
library(ggsci)
library(glue)
library(patchwork)
library(tidyverse)

font_import(paths = "C:/Windows/Fonts", pattern = "Helvetica.tff")

cr_dir <- "."
cr_file <- "cr_cytosines.csv"
hf_file <- "hf_bias_choosed_cytosines.csv"

cr_df <- read.csv(glue("{cr_dir}/{cr_file}"))
hf_df <- read.csv(glue("{cr_dir}/{hf_file}"))
print(unique(hf_df$lab))
cr_ridge_plot <- function(
        cr_dataframe,
        included_labs,
        type_filter,
        font_size = 12,
        x_limits = c(0, 10),
        x_steps = 1,
        y_label = TRUE,
        cmap = NULL
) {
    filtered_data <- cr_dataframe %>%
        filter(lab %in% included_labs,
               type == type_filter) %>%
        mutate(lab = factor(lab, levels = rev(included_labs)))
    
    p <- ggplot(filtered_data, aes(x = cr, y = lab, fill = lab))
    
    if (!is.null(cmap)) {
        if (length(cmap) != length(included_labs))
            stop("cmap 的长度必须与 included_labs 相同")
        colors_named <- setNames(rev(cmap), rev(included_labs))
        p <- ggplot(filtered_data, aes(x = cr, y = lab, fill = lab)) +
            scale_fill_manual(values = colors_named)
    } else {
        p <- ggplot(filtered_data, aes(x = cr, y = lab, fill = lab)) +
            scale_fill_npg()
    }
    
    p <- p + geom_density_ridges(
        scale = 1.5,
        alpha = 0.7,
        rel_min_height = 0.01,
        linewidth = 1.0
    ) +
        scale_x_continuous(
            limits = x_limits,
            breaks = seq(x_limits[1], x_limits[2], by = x_steps)
        ) +
        labs(title = NULL, x = NULL, y = NULL) +
        theme_ridges(
            grid = TRUE,
            center_axis_labels = TRUE
        ) +
        theme(
            text = element_text(family = "Helvetica", size = font_size),
            legend.position = "none",
            plot.title = element_blank(),
            axis.line.x = element_line(linewidth = .5, color = "black"),
            axis.line.y = element_line(linewidth = .5, color = "black"),
            axis.text.y = if (y_label) element_text(family = "Helvetica") else element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_line(linetype = "dashed", color = "gray80"),
            panel.grid.minor = element_blank(),
            legend.background =  element_blank()
        )
    return(p)
}

hf_ridge_plot <- function(
        hf_dataframe,
        included_labs,
        font_size = 12,
        x_limits = c(-60, 60),
        x_steps = 30,
        y_label = FALSE
) {
    filtered_data <- hf_dataframe %>%
        filter(lab %in% included_labs,
               direction != "no bias") %>%
        mutate(lab = factor(lab, levels = rev(included_labs)))
    
    p <- ggplot(filtered_data, aes(x = bias, y = lab, fill = direction)) +
      
        geom_density_ridges(
            scale = 1.5,
            alpha = 0.7,
            rel_min_height = 0.01,
            linewidth = 1.0
        ) + 
        geom_vline(xintercept = 0, linetype = "dashed",
                   color = "gray40", linewidth = 1.0) +
        labs(title = NULL, x = NULL, y = NULL) +
        theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
        theme(
            text = element_text(family = "Helvetica", size = font_size),
            plot.title = element_blank(),
            axis.line.x = element_line(linewidth = 0.5, color = "black"),
            axis.line.y = element_line(linewidth = .5, color = "black"),
            axis.title.y = element_blank(),
            axis.text.y = if (y_label) element_text(family = "Helvetica") else element_blank(),
            axis.ticks = element_line(linewidth = 0.5, color = "black"),
            axis.ticks.length = unit(5, "pt"),
            panel.grid.major = element_line(linetype = "dashed", color = "gray80"),
            panel.grid.minor = element_blank(),
            legend.background =  element_blank(),
            legend.position = "none") +
        scale_fill_manual(
            values = c("under-estimate" = "#3C5488FF", "over-estimate" = "#E64B35FF"),
            name = "Direction") +
        scale_x_continuous(
            limits = x_limits,
            breaks = seq(x_limits[1], x_limits[2], by = x_steps))
    return(p)
}

bs_p1 <- cr_ridge_plot(
    cr_dataframe = cr_df,
    included_labs = c("BS1", "BS2", "BS3", "BS4"),
    cmap = c("#00A087", "#007D76", "#005965",  "#003554"),
    type_filter = "conversion insufficency",
    x_limits = c(0, 5), x_steps = 2.5, y_label = TRUE)
bs_p2 <- cr_ridge_plot(
    cr_dataframe = cr_df,
    included_labs = c("BS1", "BS2", "BS3", "BS4"),
    cmap = c("#F09386", "#ED7B6B", "#E96350", "#E64B35"),
    type_filter = "excessive conversion",
    x_limits = c(0, 10), x_steps = 5, y_label = FALSE)
bs_p3 <- hf_ridge_plot(
    hf_dataframe = hf_df,
    included_labs = c("BS1", "BS2", "BS3", "BS4"),
    x_limits = c(-60, 60),
    x_steps = 30)

# bs_ridges <- bs_p1 +
#     plot_spacer() +
#     bs_p2 +
#     plot_spacer() +
#     bs_p3 + 
#     plot_layout(ncol = 5,
#                 widths = c(1, -.1, 1, -.1, 2),
#                 guides = "collect")


em_p1 <- cr_ridge_plot(
    cr_dataframe = cr_df,
    included_labs = c("EM1", "EM2", "EM3", "EM4"),
    cmap = c("#00A087", "#007D76", "#005965",  "#003554"),
    type_filter = "conversion insufficency",
    x_limits = c(0, 5), x_steps = 2.5, y_label = TRUE)
em_p2 <- cr_ridge_plot(
    cr_dataframe = cr_df,
    included_labs = c("EM1", "EM2", "EM3", "EM4"),
    cmap = c("#F09386", "#ED7B6B", "#E96350", "#E64B35"),
    type_filter = "excessive conversion",
    x_limits = c(0, 8), x_steps = 4, y_label = FALSE)
em_p3 <- hf_ridge_plot(
    hf_dataframe = hf_df,
    included_labs = c("EM1", "EM2", "EM3", "EM4"),
    x_limits = c(-40, 40),
    x_steps = 20)

# em_ridges <- em_p1 +
#     plot_spacer() +
#     em_p2 +
#     plot_spacer() +
#     em_p3 + 
#     plot_layout(ncol = 5,
#                 widths = c(1, -.1, 1, -.1, 2),
#                 guides = "collect")

ps_p1 <- cr_ridge_plot(
    cr_dataframe = cr_df,
    included_labs = c("RR1", "GM1", "GM2", "GM3"),
    cmap = c("#00A087", "#007D76", "#005965",  "#003554"),
    type_filter = "conversion insufficency",
    x_limits = c(0, 12), x_steps = 6, y_label = TRUE)
ps_p2 <- cr_ridge_plot(
    cr_dataframe = cr_df,
    included_labs = c("RR1", "GM1", "GM2", "GM3"),
    cmap = c("#F09386", "#ED7B6B", "#E96350", "#E64B35"),
    type_filter = "excessive conversion",
    x_limits = c(0, 2), x_steps = 1, y_label = FALSE)
ps_p3 <- hf_ridge_plot(
    hf_dataframe = hf_df,
    included_labs = c("RR1", "GM1", "GM2", "GM3"),
    x_limits = c(-60, 60),
    x_steps = 30)

# ps_ridges <- ps_p1 +
#     plot_spacer() +
#     ps_p2 +
#     plot_spacer() +
#     ps_p3 + 
#     plot_layout(ncol = 5,
#                 widths = c(1, -.1, 1, -.1, 2),
#                 guides = "collect")

# rr_p1 <- cr_ridge_plot(
#     cr_dataframe = cr_df,
#     included_labs = c("RR1"),
#     cmap = c("#003554"),
#     type_filter = "conversion insufficency",
#     x_limits = c(0, 1), x_steps = .5, y_label = FALSE)
# rr_p3 <- hf_ridge_plot(
#     hf_dataframe = hf_df,
#     included_labs = c("RR1"),
#     x_limits = c(-60, 60),
#     x_steps = 30)
# 
# rr_ridges <- (rr_p1 + rr_p3) + 
#     plot_layout(ncol = 2, nrow = 1,
#                 widths = c(2, 1),
#                 guides = "collect")


rm_ma_p3 <- hf_ridge_plot(
    hf_dataframe = hf_df,
    included_labs = c("RM1", "MA1", "MA2", "MA3", "NP1"),
    x_limits = c(-60, 60),
    x_steps = 30,
    y_label = TRUE)
rm_ma_p3

layout <- "
AABBCCCC#DDEEFFFF
AABBCCCC#DDEEFFFF
AABBCCCC#DDEEFFFF
AABBCCCC#DDEEFFFF
AABBCCCC#DDEEFFFF
AABBCCCC#DDEEFFFF
AABBCCCC#DDEEFFFF
AABBCCCC#DDEEFFFF
#################
GGHHIIII#####JJJJ
GGHHIIII#####JJJJ
GGHHIIII#####JJJJ
GGHHIIII#####JJJJ
GGHHIIII#####JJJJ
GGHHIIII#####JJJJ
GGHHIIII#####JJJJ
GGHHIIII#####JJJJ
"
cr_bias <- (
    bs_p1 + bs_p2 + bs_p3 +
        em_p1 + em_p2 + em_p3 +
        ps_p1 + ps_p2 + ps_p3 +
        rm_ma_p3) +
    plot_layout(
        design = layout,
        guides = "collect")

# ggsave(glue("{cr_dir}/cr_bias.png"), cr_bias,
#        width = 12,height = 6)
topptx(cr_bias, glue("{cr_dir}/cr_bias_new.pptx"),
       width = 12,height = 6)

cr_bias
