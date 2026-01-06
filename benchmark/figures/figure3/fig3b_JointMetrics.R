library(tidyverse)
library(scales)
library(linkET)
library(vegan)
library(glue)
library(eoffice)

accuracy_dir <- "."
metrics <- read.csv(glue("{accuracy_dir}/joined_metrics.csv")) %>%
  as_tibble() %>%
  column_to_rownames("sample")

# 右侧：RMSE 和 Overlap_Reference
metrics_main <- metrics %>%
  select(RMSE, Overlap_Reference) %>%
  as.data.frame()

# 左侧：其余指标
# metrics_others <- metrics %>%
#   select(-RMSE, -Overlap_Reference) %>%
#   as.data.frame()
metrics_others <- metrics %>%
    select(!all_of(c(
        "RMSE",
        "Overlap_Reference",
        "depth.5x...",
        "Q20...",
        "GC.",
        "Mean_Reads_Length",
        "Ambiguous_Characters"
    ))) %>%
    as.data.frame()

mantel <- mantel_test(metrics_main, metrics_others,
                      spec_select = list(RMSE=1,Overlap_Reference=2))%>% 
  mutate(rd = cut(r, breaks = c(-Inf,  0.5, Inf),
                  labels = c("< 0.5", ">= 0.5")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

original_names <- colnames(metrics_others)
cor_data <- correlate(metrics_others,
                      method = "spearman",
                      use = "pairwise.complete.obs")
attr(cor_data$r, "dimnames") <- list(original_names, original_names)
attr(cor_data$p, "dimnames") <- list(original_names, original_names)

q <- qcorrplot(correlate(metrics_others,
                         method = "spearman",
                         use = "pairwise.complete.obs"), 
          type = "lower", diag = FALSE) +
  geom_square() +geom_mark(sep = '\n',size = 1.8, sig_level = c(0.05, 0.01, 0.001),
                           sig_thres = 0.05,color="white") +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02"  , "#1B9E77"  , "#CCCCCC99")) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(color = "black"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3)) +
  theme(text = element_text(family = "Arial"))
q

topptx(q, glue("{accuracy_dir}/joint_metrics_new.pptx"),
       width = 10, height = 8)
