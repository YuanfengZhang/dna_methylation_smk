library(tidyverse)
library(ggforce)
library(patchwork)
library(extrafont)
library(eoffice)
library(glue)

plot_methylation_composition <- function(values, colors = c("#3C5488", "#DC0000", "#F39B7F")) {
    # 输入验证
    if (length(values) != 3) stop("values must be a numeric vector of length 3: c(5mC, 5hmC, mixed)")
    if (!is.numeric(values)) stop("values must be numeric")
    values <- pmax(values, 0)
    total <- sum(values)
    if (total == 0) stop("Sum of values is zero. Nothing to plot.")
    
    v_5mc    <- values[1]
    v_5hmC   <- values[2]
    v_mixed  <- values[3]
    
    # 颜色验证
    if (length(colors) != 3) stop("colors must be a vector of 3 HEX colors")
    
    # ========================
    # 主数据：5mC vs Others
    # ========================
    others_pct <- (v_5hmC + v_mixed) / total * 100
    others_span_rad <- others_pct / 100 * 2 * pi
    
    # Others 中心在6点方向 如果要3点方向, 则使用pi / 2
    others_center_rad <- pi
    
    df_main <- data.frame(
        group = c("5mC", "Others"),
        value = c(100 - others_pct, others_pct)
    ) %>%
        mutate(
            angle_start = ifelse(group == "Others",
                                 others_center_rad - others_span_rad / 2,
                                 others_center_rad + others_span_rad / 2),
            angle_end = ifelse(group == "Others",
                               others_center_rad + others_span_rad / 2,
                               2*pi + others_center_rad - others_span_rad / 2),
            x0 = 0, y0 = 0, r0 = 0.6, r1 = 1
        ) %>%
        mutate(group = as.character(group))  # 确保是字符型
    
    # 5mC 标签
    df_main_label <- data.frame(
        x = 0,
        y = 0.8,
        label = sprintf("%.1f%%", 100 - others_pct)  # 5mC
    )
    
    # ========================
    # 主图
    # ========================
    p_main <- ggplot() +
        geom_arc_bar(aes(x0 = x0, y0 = y0, r0 = r0, r = r1,
                         start = angle_start, end = angle_end, fill = group),
                     data = df_main) +
        geom_text(data = df_main_label,
                  aes(x = x, y = y, label = label),
                  size = 3.5, color = "white", fontface = "bold") +
        coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
        scale_fill_manual(
            values = c(colors[1], "lightgray"),  # 5mC -> colors[1], Others -> lightgray
            limits = c("5mC", "Others"),         # 强制顺序
            guide = "none"
        ) +
        theme_void() +
        theme(, plot.margin = margin(0, 0, 0, 0))
    
    # ========================
    # 小数据：5hmC 和 mixed
    # ========================
    sub_values <- c(v_5hmC, v_mixed)
    sub_total <- sum(sub_values)
    if (sub_total == 0) {
        sub_spans <- rep(pi, 2)
    } else {
        sub_spans <- sub_values / sub_total * 2 * pi
    }
    
    # 旋转偏移：+ pi/2（90度）
    rotation_offset <- pi/2
    
    df_sub <- data.frame(
        group = c("5hmC", "mixed"),
        value = sub_values
    ) %>%
        mutate(
            center_angle = ifelse(group == "mixed", pi/2, 3*pi/2) + rotation_offset,
            angle_start = center_angle - sub_spans / 2,
            angle_end = center_angle + sub_spans / 2,
            x0 = 0, y0 = 0, r0 = 0.6, r1 = 1
        ) %>%
        mutate(group = as.character(group))
    
    # 小圆环标签
    df_sub_labels <- data.frame(
        group = c("5hmC", "mixed"),
        x = 0,
        y = c(-0.8, 0.8),
        label = c(
            sprintf("%.1f%%", sub_values[2]),   # mixed
            sprintf("%.1f%%", sub_values[1])    # 5hmC
        )
    )
    
    # ========================
    # 小图
    # ========================
    p_sub <- ggplot() +
        geom_arc_bar(aes(x0 = x0, y0 = y0, r0 = r0, r = r1,
                         start = angle_start, end = angle_end, fill = group),
                     data = df_sub) +
        geom_text(data = df_sub_labels,
                  aes(x = x, y = y, label = label),
                  size = 3, color = "white", fontface = "bold") +
        coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
        scale_fill_manual(
            values = c(colors[2], colors[3]),    # 5hmC -> colors[2], mixed -> colors[3]
            limits = c("5hmC", "mixed"),
            guide = "none"
        ) +
        theme_void() +
        theme(plot.margin = margin(0, 0, 0, 0))
    
    # ========================
    # 拼接
    # ========================
    plot <- (p_main / plot_spacer() / p_sub) + 
        plot_layout(heights = c(1, -0.1, 1))
    
    return(plot)
}


setwd(".")
data_dir <- "."
output_pptx <- glue("{data_dir}/nanopore_stat.pptx")

nanopore_df <- read_csv(glue("{data_dir}/nanopore_stat.csv"),
                        col_types = cols(
                            label = col_character(),
                            total = col_integer(),
                            `annotated%` = col_double(),
                            `5mC%` = col_double(),
                            `5hmC%` = col_double(),
                            `mixed%` = col_double()
                        ))

for (i in seq_len(nrow(nanopore_df))) {
    row <- nanopore_df[i, ]
    
    # 提取数值
    values <- c(row$`5mC%`, row$`5hmC%`, row$`mixed%`)
    sample_name <- row$label
    
    # 创建图形
    p <- plot_methylation_composition(
        values = values,
        colors = c("#3498db", "#e67e22", "#c0392b")  # 可自定义
    )
    
    # 添加标题（可选）
    p <- p + ggplot2::labs(title = sample_name)
    
    # 保存到 PPTX（追加模式）
    eoffice::topptx(
        p,
        file = output_pptx,
        append = TRUE # 关键：追加到同一个文件
    )
    message("已保存：", sample_name)
}

message("✅ 所有图表已保存至：", output_pptx)