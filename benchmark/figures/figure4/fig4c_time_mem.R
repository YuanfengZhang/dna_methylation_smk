library(tidyverse)
library(eoffice)
library(readxl)
library(cowplot)

time <- read_excel("formatted_pipelines.xlsx")
# aligner
BS_ali <- time %>%
  filter(stage == "aligner") 
BS_ali_mean <- BS_ali %>%
  group_by(aligner) %>%
  summarise(time = mean(`time(s/GB)`, na.rm = T),
            rss = mean(`peak_RSS(MB/GB)`, na.rm = T)) %>%
  arrange(rss) %>%
  mutate(aligner = factor(aligner, levels = aligner))

p_time <- ggplot(BS_ali_mean, aes(x = aligner, y = time)) +
  geom_col(fill = "#1f78b4") +
  labs(x = "Time (s/GB)", y = "Aligner") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid   = element_blank() 
  )
p_time
p_rss <- ggplot(BS_ali_mean, aes(x = aligner, y = rss)) +
  geom_col(fill = "#33a02c") +
  labs(x = "RSS (GB)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10), 
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid   = element_blank() 
  )
p_rss
p_combined <- plot_grid(
  p_time, p_rss,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(1, 1)   
);p_combined


topptx(p_combined,"aligner_time.pptx",width = 8,height = 6)

# counter
BS_count <- time %>%
  filter(stage == "counter") 
BS_count_mean <- BS_count %>%
  group_by(counter) %>%
  summarise(time = mean(`time(s/GB)`, na.rm = T),
            rss = mean(`peak_RSS(MB/GB)`, na.rm = T))%>%
  arrange(rss) %>%
  mutate(counter = factor(counter, levels = counter))

p_time <- ggplot(BS_count_mean, aes(x = counter, y = time)) +
  geom_col(fill = "#1f78b4") +
  labs(x = "Time (s/GB)", y = "counter") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid   = element_blank() 
  )
p_time
p_rss <- ggplot(BS_count_mean, aes(x = counter, y = rss)) +
  geom_col(fill = "#33a02c") +
  labs(x = "RSS (GB)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10),   
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid   = element_blank() 
  )
p_rss
p_combined <- plot_grid(
  p_time, p_rss,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(1, 1)   
);p_combined


topptx(p_combined,"counter_time.pptx",width = 8,height = 6)

# deduper
BS_de <- time %>%
  filter(stage == "deduper") 
BS_de_mean <- BS_de %>%
  group_by(deduper) %>%
  summarise(time = mean(`time(s/GB)`, na.rm = T),
            rss = mean(`peak_RSS(MB/GB)`, na.rm = T))%>%
  arrange(rss) %>%
  mutate(deduper = factor(deduper, levels = deduper))

p_time <- ggplot(BS_de_mean, aes(x = deduper, y = time)) +
  geom_col(fill = "#1f78b4") +
  labs(x = "Time (s/GB)", y = "Deduper") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid   = element_blank() 
  )
p_time
p_rss <- ggplot(BS_de_mean, aes(x = deduper, y = rss)) +
  geom_col(fill = "#33a02c") +
  labs(x = "RSS (GB)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10), 
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid   = element_blank() 
  )
p_rss
p_combined <- plot_grid(
  p_time, p_rss,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(1, 1)   
);p_combined


topptx(p_combined,"deduper_time.pptx",width = 8,height = 6)

# cali
BS_ca <- time %>%
  filter(stage == "calibrator") 
BS_ca_mean <- BS_ca %>%
  group_by(calibrator) %>%
  summarise(time = mean(`time(s/GB)`, na.rm = T),
            rss = mean(`peak_RSS(MB/GB)`, na.rm = T))%>%
  arrange(rss) %>%
  mutate(calibrator = factor(calibrator, levels = calibrator))

p_time <- ggplot(BS_ca_mean, aes(x = calibrator, y = time)) +
  geom_col(fill = "#1f78b4") +
  labs(x = "Time (s/GB)", y = "Calibrator") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid   = element_blank() 
  )
p_time
p_rss <- ggplot(BS_ca_mean, aes(x = calibrator, y = rss)) +
  geom_col(fill = "#33a02c") +
  labs(x = "RSS (GB)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 10), 
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    panel.grid   = element_blank() 
  )
p_rss
p_combined <- plot_grid(
  p_time, p_rss,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(1, 1)   
);p_combined


topptx(p_combined,"calibrator_time.pptx",width = 8,height = 6)
