library(lme4)
library(splines)
library(dplyr)
library(ggplot2)
library(patchwork)
library(eoffice)


data_path <- "grouped_rmse"
df <- read.csv(data_path,
               stringsAsFactors = FALSE,
               check.names = FALSE)
processed_data <- df %>%
    filter(over_under %in% c("under", "over"),
           lab != "BS4") %>%
    mutate(
        meth_label = methylation_level,
        meth_mid = case_when(
            methylation_level == "0–10"   ~ 5,
            methylation_level == "10–20"  ~ 15,
            methylation_level == "20–30"  ~ 25,
            methylation_level == "30–40"  ~ 35,
            methylation_level == "40–50"  ~ 45,
            methylation_level == "50–60"  ~ 55,
            methylation_level == "60–70"  ~ 65,
            methylation_level == "70–80"  ~ 75,
            methylation_level == "80–90"  ~ 85,
            methylation_level == "90–100" ~ 95,
            TRUE                           ~ as.numeric(methylation_level)
        ),
        lab = factor(lab),
        sample = factor(sample),
        rep = factor(rep),
        over_under = factor(over_under, levels = c("under", "over")),
        log_count = log(count),
        method = factor(substr(lab, 1, 2))
    ) %>%
    filter(!is.na(meth_mid)) %>%
    select(lab, sample, rep, meth_mid, meth_label,
           over_under, count, log_count, method, rmse)

processed_data <- processed_data %>%
    mutate(binom_var = sqrt((meth_mid/100) * (1 - meth_mid/100)))

mixed_lm <- lmer(
    rmse ~ meth_mid * over_under +
        (1 | method) +
        (1 | lab) +
        (1 | sample),
    weights = log_count,
    data = processed_data,
    REML = TRUE
)
print(summary(mixed_lm))

mixed_lm_binom <- lmer(
  rmse ~ binom_var * over_under +
    (1 | method) +
    (1 | lab) +
    (1 | sample),
  weights = log_count,
  data = processed_data,
  REML = TRUE
)
print(summary(mixed_lm_binom))

spline_lm <- lmer(
  rmse ~ ns(meth_mid, df=3) * over_under +
    (1 | method) +
    (1 | lab) +
    (1 | sample),
  weights = log_count,
  data = processed_data,
  REML = TRUE
)
print(summary(spline_lm))

spline_lm_binom <- lmer(
  rmse ~ ns(binom_var, df=3) * over_under +
    (1 | method) +
    (1 | lab) +
    (1 | sample),
  weights = log_count,
  data = processed_data,
  REML = TRUE
)
print(summary(spline_lm_binom))

assay_linechart <- function(data, assay) {
    plot_data <- subset(data, lab == assay)

    labels_order <- plot_data %>%
        distinct(meth_label, meth_mid) %>%
        arrange(meth_mid) %>%
        pull(meth_label)

    plot_data$meth_label <- factor(plot_data$meth_label,
                                   levels = labels_order)

    p <- ggplot(data = plot_data,
                aes(x = meth_label, y = rmse,
                    color = over_under)) +
        geom_boxplot(
            position = position_dodge(width = 0.8),
            width = 0.5,
            linewidth = 0.6,
            outlier.shape = NA,
            alpha = 0.4) +
        geom_line(
            stat = "summary",
            fun = median,
            linewidth = 1,
            alpha = 0.4,
            aes(group = over_under)
        ) +
        scale_color_manual(values = c("under" = "#7788ac",
                                      "over" = "#ed8172"),
                           labels = c("under", "over")) +
        scale_x_discrete(drop = FALSE) +
        theme_bw(base_family = "Arial",
                 base_size = 12) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey90",
                                            size = 0.3),
            legend.position = "none",
            axis.ticks = element_line(size = 0.6,
                                      color = "black"),
            axis.text.x = element_text(angle = 45,
                                       hjust = 1)) +
        labs(x = "methylation level (%)",
             y = "RMSE", fill = NULL)

    return(p)
}

assay_list <- c("BS1", "BS2", "BS3", "BS4",
                "EM1", "EM2", "EM3", "EM4",
                "RR1", "GM1", "GM2", "GM3")
plot_list <- lapply(assay_list, function(assay_name) {
    assay_linechart(data = processed_data, assay = assay_name)
})
final_patch <- wrap_plots(plot_list, ncol = 2, guides = "collect",
                          axes = "collect_x",
                          axes_titles = "collect")
topptx(final_patch, "test.pptx", width = 8.5, height = 11)
