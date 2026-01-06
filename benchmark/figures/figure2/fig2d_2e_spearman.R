# 2025-0820
# EQA unified_spearman
library(gghalves)
library(ggsci)
library(ggh4x)
library(eoffice)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)

rep <- read.table("unified_spearmanr.csv",sep = ",",header = T)
head(rep)
rep1 <- rep[rep$feature=="global",]
rep1 <- rep1[,c(1,2,6)]
rep1 <- unique(rep1)
head(rep1)
print(unique(rep1$lab))
rep1_data <- rep1 %>%
  pivot_wider(names_from = label, values_from = spearmanr)
head(rep1_data)

rep1_data<-as.data.frame(rep1_data)
rownames(rep1_data)<-rep1_data$lab
rep1_data<-rep1_data[,-1]
rep1_data <- rep1_data %>%
  select(BC,BL,D5,D6,F7,M8,T1,T2,T3,T4,HF)
head(rep1_data)
# anno_20
anno_20 <- as.data.frame(rownames(rep1_data))
colnames(anno_20)<-"Filename"
anno_20$Platform <- sapply(strsplit(as.character(anno_20$Filename),""),function(x){paste(x[1],x[2],sep = "")})

anno_20$Platform <- gsub("BS","WGBS",anno_20$Platform)
anno_20$Platform <- gsub("EM","EM-seq",anno_20$Platform)
anno_20$Platform <- gsub("GM","GM-seq",anno_20$Platform)
anno_20$Platform <- gsub("RR","RRBS",anno_20$Platform)
anno_20$Platform <- gsub("RM","RRMS",anno_20$Platform)

rownames(anno_20)<-anno_20$Filename
# anno_20<-anno_20[,-1]
anno_20$Filename <- NULL
# anno_12
anno_12 <- as.data.frame(colnames(rep1_data))
colnames(anno_12)<-"Filename"
anno_12$Type <- sapply(strsplit(as.character(anno_12$Filename),""),function(x){paste(x[1],x[2],sep = "")})
rownames(anno_12)<-anno_12$Filename
# anno_12<-anno_12[,-1]
anno_12$Filename <- NULL

ann_colors3= list(
  Platform=c(WGBS = "#FFDDAD",
             "EM-seq" = "#AEC5EB",
             RRBS = "#049A8F",
             "GM-seq" = "#3A405A",
             RRMS = "#BEE5A0",
             MA= "#3C5487"), 
  Type = c(
    D5 = "#4CC3D9",
    D6 = "#7BC8A4",
    F7 = "#FFC65D",
    M8 = "#F16745",
    T1 = "gray100",
    T2 = "gray66",
    T3 = "gray33",
    T4 = "gray0",
    BC = "#69001f",
    BL = "#f7b293",
    HF ="#8039f9"
  ))  


# 按行总和从大到小排序
row_sums <- rowSums(rep1_data, na.rm = TRUE)
rep1_data_sorted <- rep1_data[order(row_sums, decreasing = TRUE), ]
rownames(rep1_data_sorted)
lab_order <- c("MA3","MA2","MA1","RR1","GM3","GM1","GM2","EM3","BS2","EM2", "EM4","BS3", "RM1","EM1","BS4","BS1")
sample_order <- c("D5","D6","F7","M8","T1","T2","T3","T4","BC","BL","HF")
rep1_data1 <- rep1_data[lab_order,sample_order]
rep_p1 <- pheatmap(rep1_data1,
                    cluster_cols = F,
                    cluster_rows = F,
                    annotation_legend = T,
                    treeheight_row = 12,treeheight_col = 12,
                    scale = "none",
                    # scale = "column",
                    # gaps_row = c(12),
                    # gaps_col = c(24, 42),
                    breaks = seq(0, 1, by = 0.1),
                    color = colorRampPalette(c("gray33","white",'#FFC65D'))(10),
                    cellwidth = 14, cellheight = 14,
                    clustering_distance_cols="correlation",
                    clustering_distance_rows="correlation",
                    #clustering_distance_cols = "euclidean",
                    #clustering_distance_rows = "euclidean",
                    show_colnames =F,show_rownames = T,
                    annotation_row=anno_20,
                    annotation_col=anno_12,
                    #  border=F,
                    border_color = "black",
                    annotation_colors = ann_colors3);rep_p1

topptx(rep_p1,"spearman.pptx",width = 8,height = 6)

min(rep1_data1)
max(rep1_data1)

# protocol_spearmanr
sp_p <- rep1 %>%
  mutate(protocol = substr(lab,1,2))
sp_p$protocol <- gsub("BS","WGBS",sp_p$protocol)
sp_p$protocol <- gsub("EM","EM-seq",sp_p$protocol)
sp_p$protocol <- gsub("GM","GM-seq",sp_p$protocol)
sp_p$protocol <- gsub("RR","RRBS",sp_p$protocol)
sp_p$protocol <- gsub("RM","RRMS",sp_p$protocol)
sp_p$protocol <- gsub("NP","Nanopore",sp_p$protocol)

sp_p$protocol <- factor(sp_p$protocol,
                          levels = c("WGBS","EM-seq","RRBS","GM-seq","RRMS","Nanopore","MA"))
sp_p$lab <- factor(sp_p$lab,
                     levels = c("BS1", "BS2", "BS3", "BS4",
                                "EM1", "EM2", "EM3", "EM4",
                                "RR1",
                                "GM1", "GM2", "GM3",
                                "RM1",
                                "MA1", "MA2", "MA3"))
sp_p$label <- factor(sp_p$label,
                   levels = c("D5","D6","F7","M8","T1","T2","T3","T4","BC","BL","HF"))
method_colors <- c(
  "D5"='#4CC3D9',"D6"='#7BC8A4',"F7"='#FFC65D',"M8"='#F16745',
  "BC"="#69001f","BL"="#f7b293",
  "T1"="gray100","T2"="gray66","T3"="gray33","T4"="gray0",
  "HF"="#8039f9",
  "BS1" = "#FFDDAD","BS2" = "#FFDDAD","BS3" = "#FFDDAD","BS4" = "#FFDDAD",
  "EM1" = "#AEC5EB","EM2" = "#AEC5EB","EM3" = "#AEC5EB","EM4" = "#AEC5EB",
  "RR1" = "#049A8F",
  "GM1" = "#3A405A","GM2" = "#3A405A","GM3" = "#3A405A",
  "RM1" = "#BEE5A0",
  "MA1"= "#3C5487","MA2"= "#3C5487","MA3"= "#3C5487"
)
protocol_colors <- c("WGBS" = "#FFDDAD",
                     "EM-seq" = "#AEC5EB",
                     "RRBS" = "#049A8F",
                     "GM-seq" = "#3A405A",
                     "RRMS" = "#BEE5A0",
                     "MA"= "#3C5487")
p <- ggplot(data =sp_p,       
            aes(x=lab, y=spearmanr, fill=lab)) +
  geom_boxplot(
    outlier.shape = 21,           
    outlier.colour = "grey90",    
    outlier.fill = NA,            
    outlier.size = 1.2,           
    outlier.stroke = 0.85,        
    width = 0.6, 
    linewidth = 0.5
  ) +
  stat_boxplot(geom = "errorbar", 
               width = 0.6,       
               linewidth = 0.5) +
  scale_fill_manual(values = method_colors)+
  theme_bw()+
  theme(
    axis.text.x = element_text(face = "bold",color = "black"),
    axis.text.y = element_text(face = "bold",color = "black"),
    axis.title  = element_text(face = "bold",color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor  =  element_blank(),
    legend.background =  element_blank(),
    legend.title = element_text(face = "bold",color = "black"),
    legend.text = element_text(color = "black")
  )+
  facet_grid2(.~protocol,scales = "free",space="free_x",
              strip  = strip_nested(
                background_x = elem_list_rect(fill = protocol_colors,color = NA),
                text_x = element_text(color ="black",face = "bold")))+
  labs(y="Spearmanr",x="Lab");p
topptx(p,"protocol_spearmanr.pptx",width = 8,height = 4)

# inter concordance

inter <- read.table("consistency.csv",sep = ",",header = T)
head(inter)

inter_sp <- inter %>%
  filter (label == "D6") %>%
  select(lab1, lab2, spearmanr) 

flip_sp <- inter_sp %>%
      rename(lab1_tmp = lab1,
             lab2_tmp = lab2) %>%
  transmute(lab1 = lab2_tmp,
            lab2 = lab1_tmp,
            spearmanr)

labs <- unique(c(inter_sp$lab1, inter_sp$lab2))
diag_sp <- data.frame(
  lab1 = labs,
  lab2 = labs,
  spearmanr = 1
)

map_sp <- bind_rows(inter_sp,flip_sp, diag_sp) %>%
  distinct(lab1,lab2, .keep_all = TRUE)

lab1_order <- c("BS1", "BS2", "BS3", "BS4",
                "EM1", "EM2", "EM3", "EM4",
                "RR1","GM1", "GM2", "GM3",
                "RM1","NP1","MA1","MA2","MA3")
map_sp$lab1 <- factor(map_sp$lab1, levels = lab1_order)
map_sp$lab2 <- factor(map_sp$lab2, levels = lab1_order)

# pheatmap spearmanr
matrix_sp <- map_sp %>%
  pivot_wider(names_from =lab2, values_from = spearmanr) %>%
  as.data.frame()
rownames(matrix_sp) <- matrix_sp$lab1
matrix_sp <- matrix_sp[,-1]
matrix_sp <- matrix_sp[lab1_order,lab1_order]

anno_16 <- as.data.frame(rownames(matrix_sp))
colnames(anno_16) <- "Filename"
anno_16$method <- sapply(strsplit(as.character(anno_16$Filename),""),function(x){paste(x[1],x[2],sep ="")})
rownames(anno_16)<-anno_16$Filename
anno_16$Filename <- NULL
anno_16$method  <- gsub("BS","WGBS",anno_16$method)
anno_16$method  <- gsub("EM","EM-seq",anno_16$method)
anno_16$method  <- gsub("GM","GM-seq",anno_16$method)
anno_16$method  <- gsub("RR","RRBS",anno_16$method)
anno_16$method  <- gsub("RM","RRMS",anno_16$method)
anno_16$method  <- gsub("MA","935K",anno_16$method)
anno_16$method  <- gsub("NP","Nanopore",anno_16$method)


ann_colors7= list(
    method=c(WGBS = "#FFDDAD",
             "EM-seq" = "#AEC5EB",
             RRBS = "#049A8F",
             "GM-seq" = "#3A405A",
             RRMS = "#BEE5A0",
             Nanopore="#A983C6",
             "935K"= "#3C5487")) 

p <- pheatmap(matrix_sp,
              cluster_cols = F,
              cluster_rows = F,
              annotation_legend =T,
              gaps_col =c(4,8,9,12,13,16),
              treeheight_row = 16, treeheight_col = 16,
              scale ="none",
              braeks = seq(0.7,1,by = 0.05),
              color = colorRampPalette(c("#089B74", "#6AC3AB", "#CDEBE3", "#FADFE7", "#F2A1BA", "#EA638C"))(6),
              cellwidth = 14, cellheight = 14,
              show_colnames =F,show_rownames = T,
              annotation_row=anno_16,
              annotation_col=anno_16,
              border_color = "black",
              annotation_colors = ann_colors7);p

topptx(p,"inter_sp.pptx",width = 6,height = 4)
max(map_sp$spearmanr)

# pheatmap jaccard

inter_ja <- inter %>%
  filter (label == "D6") %>%
  select(lab1, lab2, jaccard) 

flip_ja <- inter_ja %>%
  rename(lab1_tmp = lab1,
         lab2_tmp = lab2) %>%
  transmute(lab1 = lab2_tmp,
            lab2 = lab1_tmp,
            jaccard)

labs <- unique(c(inter_ja$lab1, inter_ja$lab2))
diag_ja <- data.frame(
  lab1 = labs,
  lab2 = labs,
  jaccard = 1
)

map_ja <- bind_rows(inter_ja,flip_ja, diag_ja) %>%
  distinct(lab1,lab2, .keep_all = TRUE)

lab1_order <- c("BS1", "BS2", "BS3", "BS4",
                "EM1", "EM2", "EM3", "EM4",
                "RR1","GM1", "GM2", "GM3",
                "RM1","NP1","MA1","MA2","MA3")
map_ja$lab1 <- factor(map_ja$lab1, levels = lab1_order)
map_ja$lab2 <- factor(map_ja$lab2, levels = lab1_order)

matrix_ja <- map_ja %>%
  pivot_wider(names_from =lab2, values_from = jaccard) %>%
  as.data.frame()
rownames(matrix_ja) <- matrix_ja$lab1
matrix_ja <- matrix_ja[,-1]
matrix_ja <- matrix_ja[lab1_order,lab1_order]

anno_16 <- as.data.frame(rownames(matrix_ja))
colnames(anno_16) <- "Filename"
anno_16$method <- sapply(strsplit(as.character(anno_16$Filename),""),function(x){paste(x[1],x[2],sep ="")})
rownames(anno_16)<-anno_16$Filename
anno_16$Filename <- NULL
anno_16$method  <- gsub("BS","WGBS",anno_16$method)
anno_16$method  <- gsub("EM","EM-seq",anno_16$method)
anno_16$method  <- gsub("GM","GM-seq",anno_16$method)
anno_16$method  <- gsub("RR","RRBS",anno_16$method)
anno_16$method  <- gsub("RM","RRMS",anno_16$method)
anno_16$method  <- gsub("MA","935K",anno_16$method)
anno_16$method  <- gsub("NP","Nanopore",anno_16$method)


ann_colors7= list(
  method=c(WGBS = "#FFDDAD",
           "EM-seq" = "#AEC5EB",
           RRBS = "#049A8F",
           "GM-seq" = "#3A405A",
           RRMS = "#BEE5A0",
           Nanopore="#A983C6",
           "935K"= "#3C5487")) 

p <- pheatmap(matrix_ja,
              cluster_cols = F,
              cluster_rows = F,
              annotation_legend =T,
              gaps_col =c(4,8,9,12,13,16),
              treeheight_row = 16, treeheight_col = 16,
              scale ="none",
              breaks = seq(0,1,by = 0.1),
              # color = colorRampPalette(c("#089B74", "#6AC3AB", "#CDEBE3", "#FADFE7", "#F2A1BA", "#EA638C"))(6),
              color = viridis(10, option = "D"),
              cellwidth = 14, cellheight = 14,
              show_colnames =F,show_rownames = T,
              annotation_row=anno_16,
              annotation_col=anno_16,
              border_color = "black",
              annotation_colors = ann_colors7);p

topptx(p,"inter_ja.pptx",width = 8,height = 6)
max(inter_ja$jaccard)
min(map_ja$jaccard)
