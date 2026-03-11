# Load necessary packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(stats)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(glmnet)
library(caret)
library(ROCR)
library(tidyverse)
library(Boruta)
library(tidyr)
library(pheatmap)
library(gghalves)
library(RColorBrewer)
library(rstatix)
library(ggpubr)
library(circlize)
library(nnet)

# data
meta <- read.table(file = "data/all_metadata.tsv", sep = "\t", row.names = 1,
                   header = T)


data <- read.table(file = "data/605_data_unique_genenames.tsv", 
                   sep ="\t", header = T, row.names = 1)

data.pca <-merge(data, meta, by = 0)
row.names(data.pca) <- data.pca$Row.names
data.pca$Row.names <- NULL
data.pca <- data.pca[, c(1:605,626)]

# Simulating the grouping
groups <- data.pca$group # Extend this according to the actual data

data.pca$group <- factor(groups, levels = c('HC', 'MPO_R', 'PR3_R', 'MPO_A', 'PR3_A'))
# Convert to matrix
data_matrix <- as.matrix(data.pca[, -ncol(data.pca)])

# Create a data frame for row annotations
annotation_row <- data.frame(Group = factor(groups))
rownames(annotation_row) <- rownames(data_matrix)

# Define annotation colors
ann_colors <- list(
  Group = c("HC" = "darkgrey", "MPO_A" = "#cf0dcf", "MPO_R" =  "#f0b2f0",
            "PR3_A" = "#009900", "PR3_R" = "#bfe6bf")
)
# Create the heatmap
pheatmap(data_matrix,
         col = c("blue", "white", "red"), 
         cluster_rows = F, 
         cluster_cols = TRUE, 
         #clustering_distance_cols = 'euclidean',
         #clustering_distance_rows = 'euclidean',
         #clustering_method = 'ward.D',
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         annotation_names_row = FALSE, 
         annotation_names_col = FALSE,
         fontsize_row = 10,          
         fontsize_col = 7,          
         angle_col = 45, 
         legend_breaks = c(15, 30, 45), 
         legend_labels = c("Low", "Medium", "High"), 
         show_colnames = FALSE, 
         show_rownames = FALSE, 
         main = "Protein Distribution Heatmap with Clustering")



library(reshape2)
library(ggplot2)
library(RColorBrewer)

data_scaled <- scale(data_matrix)
# Define a more contrasting color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks that map the values directly from -3 to 3
breaks_list <- seq(-3, 3, length.out = 101)

# Create the heatmap with color limits between -3 and 3
pheatmap(data_scaled,
         color = color_palette, 
         breaks = breaks_list,
         scale = "none",  # Ensure the scaling is not applied since you're setting explicit limits
         cluster_rows = FALSE,  
         cluster_cols = TRUE,  
         annotation_row = annotation_row,  
         annotation_colors = ann_colors,  
         show_colnames = FALSE,
         show_rownames = FALSE,
         main = "Protein Expression Heatmap with Limits (-3 to 3)")


# Spearman correlation of median intensities of all proteins over the different groups 

data <- read.table(file = "data/605_data_unique_genenames.tsv", 
                   sep ="\t", header = TRUE, row.names = 1)

# Merge with metadata
data.pca <- merge(data, meta, by = 0)
row.names(data.pca) <- data.pca$Row.names
data.pca$Row.names <- NULL
data.pca <- data.pca[, c(1:605, 626)]

# Define groups
data.pca$group <- factor(data.pca$group, levels = c('HC', 'MPO_R', 'PR3_R', 'MPO_A', 'PR3_A'))

median_intensities <- aggregate(. ~ group, data = data.pca, median)
row.names(median_intensities) <- median_intensities$group
median_intensities$group <- NULL

spearman_corr <- cor(t(median_intensities), method = "spearman")

# Mapping original group names to custom names
group_mapping <- c('HC' = 'Healthy controls',
                   'PR3_A' = 'PR3 active',
                   'PR3_R' = 'PR3 remission',
                   'MPO_A' = 'MPO active',
                   'MPO_R' = 'MPO remission')


pheatmap::pheatmap(spearman_corr,
                   color = colorRampPalette(c("white", "red"))(100),  # Gradient from white to red
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE,
                   labels_row = group_mapping[row.names(spearman_corr)],
                   labels_col = group_mapping[colnames(spearman_corr)],
                   angle_col = 45,  # Rotate column labels by 45 degrees
                   fontsize_number = 15,  # Increase the font size of the numbers
                   fontsize = 15  # Increase the font size of the axis labels
                   # Row labels rotation is not directly supported in pheatmap
) -> pheatspearman

ggsave(pheatspearman, file = "../Uwes_wishes/April_Figure/2b_spearmanCorrelation_corrected.svg", width = 5, height = 4)

# spearmann is ranked, we dont need this here

pearson_corr <- cor(t(median_intensities), method = "pearson")

# Mapping original group names to custom names
group_mapping <- c('HC' = 'Healthy controls',
                   'PR3_A' = 'PR3 active',
                   'PR3_R' = 'PR3 remission',
                   'MPO_A' = 'MPO active',
                   'MPO_R' = 'MPO remission')

pheatmap::pheatmap(pearson_corr,
                   color = colorRampPalette(c("white", "red"))(100),  # Gradient from white to red
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE,
                   labels_row = group_mapping[row.names(pearson_corr)],
                   labels_col = group_mapping[colnames(pearson_corr)],
                   angle_col = 45,  # Rotate column labels by 45 degrees
                   fontsize_number = 15,  # Increase the font size of the numbers
                   fontsize = 15  # Increase the font size of the axis labels
                   # Row labels rotation is not directly supported in pheatmap
) -> pheatpearson_corr

ggsave(pheatpearson_corr, file = "../Uwes_wishes/April_Figure/2b_personCorrelation_corrected.svg", width = 5, height = 4)


### PCA 
## these are the 22 after the lasso regression
lasso_para <- read.table(file = "data/21_lasso_para_withgene_names.csv", sep = ",", header =T)

data <- read.table(file = "data/605_data_unique_genenames.tsv", 
                   sep ="\t", header = T, row.names = 1)

data21.pca <- data[colnames(data) %in% lasso_para$Genes]

data21.pca <-merge(data21.pca, meta, by = 0)
row.names(data21.pca) <- data21.pca$Row.names
data21.pca$Row.names <- NULL
data21.pca <- data21.pca[, c(1:20,41)]

# Dimension reduction using PCA
res.pca <- prcomp(data21.pca[, -21])


# Contributions of variables to PC1
# var = its the immunecells
factoextra::fviz_contrib(res.pca, choice = "var", axes = 1, 
             top = 10, color = "black", 
             fill = alpha( "grey", 0.5)) -> x
x + theme(#axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 90, size = 13),
  axis.text.y = element_text (size = 13), 
  axis.title.y = element_text (size = 13))->xy 

# Contributions of variables to PC2
# var = its the immunecells
factoextra::fviz_contrib(res.pca, choice = "var", axes = 2, 
             top = 10, color = "black", 
             fill = alpha( "grey", 0.5)) -> y
y + theme(#axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 90, size = 13),
  axis.text.y = element_text (size = 13), 
  axis.title.y = element_text (size = 13))->yx 
cowplot::plot_grid(xy,yx, nrow = 2) -> contribution

library(factoextra)

ind.coord.2 <- as.data.frame(get_pca_ind(res.pca)$coord)
groups <- data21.pca[,c(1,21)]
groups$LRG1 <- NULL

ind.coord <- merge(ind.coord.2, groups, by = 0)
row.names(ind.coord) <- ind.coord$Row.names
ind.coord$Row.names <- NULL
# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

## for the different diseas.groups
ind.coord.HC <- ind.coord[ind.coord$group %in% c("HC"), ]
ind.coord.MPO_A <- ind.coord[ind.coord$group %in% c("MPO_A"), ]
ind.coord.MPO_R <- ind.coord[ind.coord$group %in% c("MPO_R"), ]
ind.coord.PR3_A <- ind.coord[ind.coord$group %in% c("PR3_A"), ]
ind.coord.PR3_R <- ind.coord[ind.coord$group %in% c("PR3_R"), ]

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

col_vector <- c("darkgrey", "#cf0dcf",  "#f0b2f0","#009900","#bfe6bf")

centroids <- aggregate(cbind(Dim.1,Dim.2)~ group,data = ind.coord,mean)

ggplot(ind.coord, aes(x = Dim.1, y = Dim.2, fill = as.factor(group))) +
  theme_classic() + 
  scale_fill_manual(values = col_vector) +
  scale_color_manual(values = col_vector) +
  geom_point(aes(fill = as.factor(group)), size = 4, pch = 21, alpha = 0.8) + 
  xlab("Dim 1 (35.1%)") + ylab("Dim 2 ( 16.3%)") + 
  theme(axis.title.x = element_text(size = 13), 
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        axis.title.y = element_text(size = 13),
        legend.position = "right") +
  stat_ellipse(aes(color = as.factor(group), fill = as.factor(group)), 
               geom = "polygon", 
               alpha = 0.2) + 
  geom_point(data = centroids,
             size = 6,
             shape = 16,
             color = "black") +
  geom_point(data = centroids, 
             size = 5,pch = 21,
             shape = 16) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12)) -> a


grobs <- ggplotGrob(a)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

xdensity <- ggplot(ind.coord, aes(Dim.1)) +
  geom_density(alpha=.5, aes(color = as.factor(group), fill = as.factor(group))) +
  scale_fill_manual(values=col_vector) +
  scale_color_manual(values=col_vector) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        #axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none")

ydensity <- ggplot(ind.coord, aes(Dim.2)) +
  geom_density(alpha=.5, aes(color = as.factor(group), fill = as.factor(group))) +
  scale_fill_manual(values=col_vector) +
  scale_color_manual(values=col_vector) +
  theme_classic() +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) +
  theme(legend.position = "none") +
  coord_flip()


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )


bigpcoAHO <-
  cowplot::plot_grid(
    xdensity + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
    blankPlot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
    a + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
    ydensity + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
    nrow = 2,
    rel_widths = c(4, 1.4),
    rel_heights = c(1.4, 4),
    align = "hv"
  ) 

v.all <- cowplot::plot_grid(bigpcoAHO, legend, rel_widths = c(1, .2))

ggsave(v.all, file = "../Uwes_wishes/April_Figure/3c_21_not.scaled.pca_analysis_r.svg", 
       device = "svg", width = 6,
       height = 5)


vegan::adonis(data21.pca[, -c(21)] ~ as.factor(group), data = meta, method='canberra') -> perm.ey #0.026 *
perm.ey$aov.tab

# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# as.factor(group)   4  0.033843 0.0084608  13.524 0.25505  0.001 ***
#   Residuals        158  0.098850 0.0006256         0.74495           
# Total            162  0.132693                   1.00000           
# ---



library(factoextra)
meta <- read.table(file = "data/all_metadata.tsv", sep = "\t", row.names = 1,
                   header = T)


data <- read.table(file = "data/605_data_unique_genenames.tsv", 
                   sep ="\t", header = T, row.names = 1)

data.pca <-merge(data, meta, by = 0)
row.names(data.pca) <- data.pca$Row.names
data.pca$Row.names <- NULL
data.pca <- data.pca[, c(1:605,626)]

# Dimension reduction using PCA
library(stats)
res.pca <- prcomp(data.pca[, -606])
library(FactoMineR)
library(factoextra)

# Contributions of variables to PC1
# var = its the immunecells
fviz_contrib(res.pca, choice = "var", axes = 1, 
             top = 10, color = "black", 
             fill = alpha( "grey", 0.5)) -> x
x + theme(#axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 90, size = 13),
  axis.text.y = element_text (size = 13), 
  axis.title.y = element_text (size = 13))->xy 

# Contributions of variables to PC2
# var = its the immunecells
fviz_contrib(res.pca, choice = "var", axes = 2, 
             top = 10, color = "black", 
             fill = alpha( "grey", 0.5)) -> y
y + theme(#axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 90, size = 13),
  axis.text.y = element_text (size = 13), 
  axis.title.y = element_text (size = 13))->yx 
cowplot::plot_grid(xy,yx, nrow = 2) -> contribution


ind.coord.2 <- as.data.frame(get_pca_ind(res.pca)$coord)
groups <- meta[, c(1,21)]

ind.coord <- merge(ind.coord.2, groups, by = 0)
row.names(ind.coord) <- ind.coord$Row.names
ind.coord$Row.names <- NULL
# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

## for the different diseas.groups
ind.coord.HC <- ind.coord[ind.coord$group %in% c("HC"), ]
ind.coord.MPO_A <- ind.coord[ind.coord$group %in% c("MPO_A"), ]
ind.coord.MPO_R <- ind.coord[ind.coord$group %in% c("MPO_R"), ]
ind.coord.PR3_A <- ind.coord[ind.coord$group %in% c("PR3_A"), ]
ind.coord.PR3_R <- ind.coord[ind.coord$group %in% c("PR3_R"), ]

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

col_vector <- c("darkgrey", "#cf0dcf",  "#f0b2f0","#009900","#bfe6bf")


centroids <- aggregate(cbind(Dim.1,Dim.2)~ group,data = ind.coord,mean)

ggplot(ind.coord, aes (x = Dim.1, y = Dim.2, color = as.factor(group) )) +
  theme_classic () + 
  scale_color_manual(values = col_vector) +
  scale_fill_manual(values = col_vector) +
  geom_point (aes (color = as.factor(group)), size = 4, alpha = 0.8) + 
  xlab ("Dim 1(10.6%)") + ylab ("Dim 2(4.3%)")+ 
  theme (axis.title.x = element_text (size = 13), 
         axis.text.x = element_text (size = 13), 
         axis.text.y = element_text (size = 13), 
         axis.title.y = element_text (size = 13),
         legend.position="right") +
  stat_ellipse(aes(color = as.factor(group), fill = as.factor(group)), 
               geom = "polygon", 
               alpha = 0.2) +
  geom_point(data = centroids,
             size = 5,
             shape = 16,
             color = "black") +# centroides hinzufügen
  geom_point(data = centroids, 
             size = 4, 
             shape = 16)+
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15)) -> a

grobs <- ggplotGrob(a)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

xdensity <- ggplot(ind.coord, aes(Dim.1)) +
  geom_density(alpha=.5, aes(color = as.factor(group), fill = as.factor(group))) +
  scale_fill_manual(values=col_vector) +
  scale_color_manual(values=col_vector) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        #axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position="none")

ydensity <- ggplot(ind.coord, aes(Dim.2)) +
  geom_density(alpha=.5, aes(color = as.factor(group), fill = as.factor(group))) +
  scale_fill_manual(values=col_vector) +
  scale_color_manual(values=col_vector) +
  theme_classic() +
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()) +
  theme(legend.position = "none") +
  coord_flip()


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )


bigpcoAHO <-
  cowplot::plot_grid(
    xdensity + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
    blankPlot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
    a + theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "cm")),
    ydensity + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
    nrow = 2,
    rel_widths = c(4, 1.4),
    rel_heights = c(1.4, 4),
    align = "hv"
  ) 

v.all <- cowplot::plot_grid(bigpcoAHO, legend, rel_widths = c(1, .2))

ggsave(v.all, file = "../Uwes_wishes/April_Figure/1a_605not.scaled.pca_analysis_r.svg", device = "svg", width = 6,
       height = 4)

vegan::adonis(data ~ as.factor(group), data = meta, method='canberra') -> perm.ey #0.026 *
perm.ey$aov.tab
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# as.factor(group)   4   0.02186 0.0054651  6.8701 0.14816  0.001 ***
#   Residuals        158   0.12569 0.0007955         0.85184           
# Total            162   0.14755                   1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



### Vulcano plot using Log2Folds form Marie
# -Achse: log2 Fold Change Active vs Remission
# Y-Achse: -log10 p-value
# Welsh T-Test FDR corrected > FDR5

Marie_log2_welchT <- read.csv2("data/Marie_log2_welchT.csv")
# Create the volcano plot
# Ernrichment - inflammatory genes 
# take the data from the Genere enrichmet 
Inflammatory_markers_vomGenSetENrich <- read.csv("data/New_Inflammatory_markers_vomGenSetENrich.tsv", sep="\t")

# List of inflammation-related keywords
names <- Inflammatory_markers_vomGenSetENrich$geneID
names_vector <- unlist(strsplit(as.character(Inflammatory_markers_vomGenSetENrich$geneID), "/"))
names_vector <- trimws(names_vector)  # Remove any leading/trailing whitespace

# Assuming 'names_vector' contains your target inflammatory genes and is correctly formatted
Marie_log2_welchT <- Marie_log2_welchT %>%
  mutate(
    Inflammation = ifelse(PG.Genes %in% names_vector, "Inflammatory", "Non-inflammatory")
  )


Marie_log2_welchTfdr <- Marie_log2_welchT
row.names(Marie_log2_welchTfdr) <- make.names(Marie_log2_welchTfdr$PG.Genes, unique = T)

volcano_plot <- ggplot(Marie_log2_welchTfdr, aes(
  x = log2.FC.Active_Remission,
  y = p.value.Active_Remission._log10
)) +
  scale_color_manual(values = c("significant" = "black", "non_significant" = "gray")) +
  labs(
    x = "Log2 Fold Change active vs remission",
    y = "-log10(p-value)",
    color = "Significance",
    fill = "Protein Type"
  )+
  geom_point(
    data = subset(Marie_log2_welchTfdr, Inflammation == "Non-inflammatory"),
    aes(color = Welch.s.T.test.Significant.Active_Remission_FDR5),
    size = 3,  # Slightly smaller size for Non-inflammatory
    shape = 21,  # Using circles for both, can change if needed
    # color = "black",  # Fill dark yellow for non-inflammatory
    stroke = 1, alpha = 0.8
  ) +
    geom_point(
      data = subset(Marie_log2_welchTfdr, Inflammation == "Inflammatory"),
      aes(color = Welch.s.T.test.Significant.Active_Remission_FDR5),
      size = 3,  # Bigger size for Inflammatory
      shape = 21,  # Using circles for both, can change if needed
      fill = "#ff4700",  # Fill red for inflammatory
      stroke = 1, alpha = 0.8
    ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )

print(volcano_plot)

ggsave(volcano_plot, file = "Uwes_wishes/April_Figure/all_volcanos.svg", width = 4, height = 5)


volcano_plotmpo <- ggplot(Marie_log2_welchT, aes(
  x = log2.FC.MPO_active_Remission,
  y = p.value.MPO_active_Remission_log10
)) +
  scale_color_manual(values = c("significant" = "black", "non_significant" = "gray")) +
  labs(
    x = "Log2 Fold Change active vs remission",
    y = "-log10(p-value)",
    color = "Significance",
    fill = "Protein Type"
  )+
  geom_point(
    data = subset(Marie_log2_welchTfdr, Inflammation == "Non-inflammatory"),
    aes(color = Welch.s.T.test.Significant.MPO_active_MPO_Remission_FDR5),
    size = 3,  # Slightly smaller size for Non-inflammatory
    shape = 21,  # Using circles for both, can change if needed
    # color = "black",  # Fill dark yellow for non-inflammatory
    stroke = 1, alpha = 0.8
  ) +
  geom_point(
    data = subset(Marie_log2_welchTfdr, Inflammation == "Inflammatory"),
    aes(color = Welch.s.T.test.Significant.MPO_active_MPO_Remission_FDR5),
    size = 3,  # Bigger size for Inflammatory
    shape = 21,  # Using circles for both, can change if needed
    fill = "#ff4700",  # Fill red for inflammatory
    stroke = 1, alpha = 0.8
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )


ggsave(volcano_plotmpo, file = "Uwes_wishes/April_Figure/all_volcanosmpo.svg", width = 4, height = 5)

# Welch.s.T.test.Significant.PR3_active_PR3_Remission_FDR5

volcano_plotpr3 <- ggplot(Marie_log2_welchT, aes(
  x = log2.FC.PR3_active_Remission,
  y = p.value.PR3_active_Remission_log10
)) +
  scale_color_manual(values = c("significant" = "black", "non_significant" = "gray")) +
  labs(
    x = "Log2 Fold Change active vs remission",
    y = "-log10(p-value)",
    color = "Significance",
    fill = "Protein Type"
  )+
  geom_point(
    data = subset(Marie_log2_welchTfdr, Inflammation == "Non-inflammatory"),
    aes(color = Welch.s.T.test.Significant.PR3_active_PR3_Remission_FDR5),
    size = 3,  # Slightly smaller size for Non-inflammatory
    shape = 21,  # Using circles for both, can change if needed
    # color = "black",  # Fill dark yellow for non-inflammatory
    stroke = 1, alpha = 0.8
  ) +
  geom_point(
    data = subset(Marie_log2_welchTfdr, Inflammation == "Inflammatory"),
    aes(color = Welch.s.T.test.Significant.PR3_active_PR3_Remission_FDR5),
    size = 3,  # Bigger size for Inflammatory
    shape = 21,  # Using circles for both, can change if needed
    fill = "#ff4700",  # Fill red for inflammatory
    stroke = 1, alpha = 0.8
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )


# ggsave(volcano_plotpr3, file = "Uwes_wishes/April_Figure/all_volcanospr3.svg", width = 4, height = 5)



# 327 deps
names_327 <- read.csv("data/327_DEPs_MetadeconfR.csv", sep = "")
row.names(names_327) <- names_327$protein


# reduce the dataset
setdiff(row.names(names_327), row.names(Marie_log2_welchTfdr))
row.names(names_327) <- gsub("_", "\\.", row.names(names_327))
row.names(names_327) <- gsub("_", "\\.", row.names(names_327))
row.names(Marie_log2_welchTfdr) <- gsub("_", "\\.", row.names(Marie_log2_welchTfdr))
row.names(Marie_log2_welchTfdr) <- gsub("_", "\\.", row.names(Marie_log2_welchTfdr))

setdiff(row.names(names_327), row.names(Marie_log2_welchTfdr))

df.red327 <- Marie_log2_welchTfdr[row.names(Marie_log2_welchTfdr) %in% row.names(names_327), ]

# Welch.s.T.test.Significant.Active_Remission_FDR5

# Create the volcano plot
volcano_plot327 <- ggplot(df.red327, aes(
  x = log2.FC.Active_Remission,
  y = p.value.Active_Remission._log10
)) +
  scale_color_manual(values = c("significant" = "black", "non_significant" = "gray")) +
  labs(
    x = "Log2 Fold Change active vs remission",
    y = "-log10(p-value)",
    color = "Significance",
    fill = "Protein Type"
  )+
  geom_point(
    data = subset(df.red327, Inflammation == "Non-inflammatory"),
    aes(color = Welch.s.T.test.Significant.Active_Remission_FDR5),
    size = 3,  # Slightly smaller size for Non-inflammatory
    shape = 21,  # Using circles for both, can change if needed
    # color = "black",  # Fill dark yellow for non-inflammatory
    stroke = 1, alpha = 0.8
  ) +
  geom_point(
    data = subset(df.red327, Inflammation == "Inflammatory"),
    aes(color = Welch.s.T.test.Significant.Active_Remission_FDR5),
    size = 3,  # Bigger size for Inflammatory
    shape = 21,  # Using circles for both, can change if needed
    fill = "#ff4700",  # Fill red for inflammatory
    stroke = 1, alpha = 0.8
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )

print(volcano_plot327)

# ggsave(volcano_plot300, file = "Uwes_wishes/April_Figure/New2C_327DEP_volcanos.svg", width = 4, height = 5)


# Welch.s.T.test.Significant.MPO_active_MPO_Remission_FDR5

volcano_plotmpo327 <- ggplot(df.red327, aes(
  x = log2.FC.MPO_active_Remission,
  y = p.value.MPO_active_Remission_log10
))  +
  scale_color_manual(values = c("significant" = "black", "non_significant" = "gray")) +
  labs(
    x = "Log2 Fold Change active vs remission",
    y = "-log10(p-value)",
    color = "Significance",
    fill = "Protein Type"
  )+
  geom_point(
    data = subset(df.red327, Inflammation == "Non-inflammatory"),
    aes(color = Welch.s.T.test.Significant.MPO_active_MPO_Remission_FDR5),
    size = 3,  # Slightly smaller size for Non-inflammatory
    shape = 21,  # Using circles for both, can change if needed
    # color = "black",  # Fill dark yellow for non-inflammatory
    stroke = 1, alpha = 0.8
  ) +
  geom_point(
    data = subset(df.red327, Inflammation == "Inflammatory"),
    aes(color = Welch.s.T.test.Significant.MPO_active_MPO_Remission_FDR5),
    size = 3,  # Bigger size for Inflammatory
    shape = 21,  # Using circles for both, can change if needed
    fill = "#ff4700",  # Fill red for inflammatory
    stroke = 1, alpha = 0.8
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )

# ggsave(volcano_plotmpo327, file = "../Uwes_wishes/April_Figure/2C_327DEP_volcanosmpo.svg", width = 4, height = 5)

# Welch.s.T.test.Significant.PR3_active_PR3_Remission_FDR5

volcano_plotpr327 <- ggplot(df.red327, aes(
  x = log2.FC.PR3_active_Remission,
  y = p.value.PR3_active_Remission_log10
)) +
  scale_color_manual(values = c("significant" = "black", "non_significant" = "gray")) +
  labs(
    x = "Log2 Fold Change active vs remission",
    y = "-log10(p-value)",
    color = "Significance",
    fill = "Protein Type"
  )+
  geom_point(
    data = subset(df.red327, Inflammation == "Non-inflammatory"),
    aes(color = Welch.s.T.test.Significant.PR3_active_PR3_Remission_FDR5),
    size = 3,  # Slightly smaller size for Non-inflammatory
    shape = 21,  # Using circles for both, can change if needed
    # color = "black",  # Fill dark yellow for non-inflammatory
    stroke = 1, alpha = 0.8
  ) +
  geom_point(
    data = subset(df.red327, Inflammation == "Inflammatory"),
    aes(color = Welch.s.T.test.Significant.PR3_active_PR3_Remission_FDR5),
    size = 3,  # Bigger size for Inflammatory
    shape = 21,  # Using circles for both, can change if needed
    fill = "#ff4700",  # Fill red for inflammatory
    stroke = 1, alpha = 0.8
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )

# ggsave(volcano_plotpr327, file = "../Uwes_wishes/April_Figure/2C_327DEP_volcanospr3.svg", width = 4, height = 5)

long_df <- Inflammatory_markers_vomGenSetENrich %>%
  mutate(Term = Description) %>%
  select(Term, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  mutate(geneID = trimws(geneID)) %>%
  distinct()

binary_matrix <- table(long_df$Term, long_df$Count)

# Heatmap
pheatmap(binary_matrix,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Proteins in Inflammation-Related GO Terms",
         fontsize_row = 8, fontsize_col = 6)


# dataset prep:long format with term and protein
long_df <- Inflammatory_markers_vomGenSetENrich %>%
  mutate(Term = Description) %>%
  dplyr::select(Term, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  mutate(geneID = trimws(geneID)) %>%
  distinct()

# Plot
chordDiagram(long_df, transparency = 0.3)

#list of proteins that need to be included
must_include <- c("CRP", "S100A8", "ORM1", "TNC", "TIMP1")

# only inflammatory proteine and significant
df.inflam <- df.red300 %>%
  filter(
    Inflammation == "Inflammatory",
    Welch.s.T.test.Significant.Active_Remission_FDR5 == "significant"  
  )

df.inflam <- data.frame(df.inflam)


# sort by absolute Log2-FC
df.sorted <- df.inflam %>%
  arrange(desc(abs(log2.FC.Active_Remission)))

# select Top-N (25)
topN <- 25
df.topN <- df.sorted[1:topN, ]
df.topN$Protein <- row.names(df.topN)
df.inflam$Protein <- row.names(df.inflam)

# add 'must include' proteins
additional_rows <- df.inflam[df.inflam$Protein %in% must_include & !(df.inflam$Protein %in% df.topN$Protein), ]
df.top_combined <- bind_rows(df.topN, additional_rows)

# Filter long dataframe
df_long_20 <- long_df[long_df$geneID %in% df.top_combined$Protein, ]

# define colors
# Create a 15-color version of Set3
custom_set3_15 <- colorRampPalette(brewer.pal(12, "Set3"))(15)

# name them to match your terms
term_names <- unique(df_long_20$Term)
term_colors <- setNames(custom_set3_15[1:length(term_names)], term_names)

protein_colors <- rep("#4575B4", length(unique(df_long_20$geneID)))
names(protein_colors) <- unique(df_long_20$geneID)
grid_colors <- c(term_colors, protein_colors)

# Circos Plot
circos.clear()
chordDiagram(df_long_20,
             grid.col = grid_colors,
             transparency = 0.3,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.08))

# 9. Labels plotten
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + .1, sector.name,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex = 1, col = "black", font = 1)
}, bg.border = NA)


# Daten einlesen
df <- read.csv("data/Inflammatory_markers_vomGenSetENrich.csv", sep = ";")


status <- read.csv2("327_DEPs_MetadeconfR.csv")
status <- data.frame(status)
status$names <- row.names(status)
status$names <- gsub("\\_", "\\.", status$names)

pro_139_sig<- read.csv2("data/135_confoundR_cleaned.csv")
row.names(pro_139_sig) <- gsub("_", "\\.", row.names(pro_139_sig))

df.red_139 <- Marie_log2_welchTfdr[row.names(Marie_log2_welchTfdr) %in% row.names(pro_139_sig), ]
status$names <- gsub("_", "\\.", status$names)
status <- status[status$names %in% row.names(df.red_139), ]
row.names(status) <- status$names
comb <- merge(df.red_139, status, by = 0)
row.names(comb) <- comb$Row.names
comb$Row.names <- NULL

vulcano_input <- comb
vulcano_input$active_vs_remission_status_log10 <-log10(vulcano_input$active_vs_remission_status)
vulcano_input$active_vs_remission_status_log10  <- abs(vulcano_input$active_vs_remission_status_log10 )

parameter.after.lasso.1 <- read.csv("data/parameter.after.lasso.1.csv", header=T)

features_to_label <- parameter.after.lasso.1[c(2:22), ]
features_to_label <- unlist(strsplit(as.character(features_to_label$names), ","))
features_to_label <- trimws(features_to_label)  # Remove any leading/trailing whitespace

# Add a column to the dataframe that indicates whether the feature should be labeled
vulcano_input$label <- ifelse(row.names(vulcano_input) %in% features_to_label, vulcano_input$PG.Genes, "")
vulcano_input$sig[(vulcano_input$active_vs_remission_status <= 0.1)] <- "sig"
vulcano_input <- vulcano_input %>% 
  filter(!sig == "NA")

library(ggrepel)
volcano_plot135 <- ggplot(vulcano_input, aes(
  x = log2.FC.Active_Remission,
  y = active_vs_remission_status_log10
)) +
  geom_point(
    aes(fill = label),  # Color the dots based on the label column
    shape = 21, color = "black", alpha = 0.7, size = 3, stroke = 1
  ) + 
  scale_fill_manual(values = c(label = "darkgreen")) +
  scale_color_manual(values = c("significant" = "black", "NA" = "black")) +
  labs(
    x = "Log2 Fold Change active vs remission (135)",
    y = "-log10(p-value)",
    color = "Features"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  geom_label_repel(
    aes(label = label), 
    size = 4,
    box.padding = 0.5,
    max.overlaps = 100,
    segment.color = "black",  # Color of the segment lines to the labels
    fill = "white",  # Background color of the text label
    alpha = 0.8,  # Transparency level of the background color
    family = "arial"  # Font family, can be adjusted as needed
  )



# helper: TRUE if label is not empty
vulcano_input$has_label <- ifelse(vulcano_input$label != "", "yes", "no")

# Plot
volcano_plot135 <- ggplot(vulcano_input, aes(
  x = log2.FC.Active_Remission,
  y = active_vs_remission_status_log10
)) +
  geom_point(
    aes(fill = has_label),
    shape = 21, color = "black", alpha = 0.7, size = 4, stroke = 1
  ) +
  scale_fill_manual(values = c("yes" = "darkgreen", "no" = "white")) +  # Nur Punkte mit Label einfärben
  labs(
    x = "Log2 Fold Change active vs remission (135)",
    y = "-log10(p-value)",
    fill = "Labeled"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  geom_label_repel(
    aes(label = label),
    size = 4,
    box.padding = 0.5,
    max.overlaps = 100,
    segment.color = "black",
    fill = "white",
    alpha = 0.8,
    family = "arial"
  )

# ggsave(volcano_plot139, file = "../Uwes_wishes/April_Figure/3bNEw_volcano_plot135_with21Lasso.svg", width = 3.5, height = 4)


# AUC 21 Panel on the 20 % test of the Berlin cohort 
data <- read.table(file = "data/20test.data.csv", 
                   sep =",", header = T, row.names = 1)

# Define the mapping of old names to new names
rename_mapping <- c(
  "PR3_active_A20.40" =  "MPO_active_A20.40",
  "PR3_active_A21.14" = "MPO_active_A21.14" ,
  "PR3_active_A21.6" = "MPO_active_A21.6" ,
  "PR3_active_A20.21.1" =  "PR3_active_A20.21_a",
  "Prag_PR3_Remission_5631_P02_A02" = "Prag_PR3_Active_5631_P02_A02"
)

row.names(data) <- gsub("MPA", "MPO", row.names(data))

# Rename the row names in `data`
row.names(data) <- ifelse(row.names(data) %in% names(rename_mapping),
                          rename_mapping[row.names(data)],
                          row.names(data))



# lets try the 21 panel 
# Load metadata
meta <- read.table(file = "data/all_metadata.tsv", sep = "\t", header = T, row.names = 1)
# Data preprocessing
meta <- meta %>%
  filter(!group=="HC")
meta$Pseudonym <- NULL
meta$Sample.Number <- NULL
meta <- meta %>%
  mutate(
    sex = if_else(Sex == 'female', 0, 1)
  )
# Additional data preprocessing
row.names(meta) <- meta$Sample_ID_Proteomics
meta$Sample_ID_Proteomics <- NULL
meta$group <- NULL
meta$ANCA_specifity_Disease <- NULL
meta$disease.nr <- NULL


prepare_data <- function(data_subset) {
  data_subset$Experimental_Group <- row.names(data_subset)
  data_subset <- data_subset %>%
    mutate(group = case_when(
      grepl("active", Experimental_Group, ignore.case = TRUE) ~ 0,  # Active = 0
      grepl("remission", Experimental_Group, ignore.case = TRUE) ~ 1,  # Remission = 1
      TRUE ~ NA_real_
    ))
  data_subset$group <- as.factor(data_subset$group)
  return(data_subset)
}

meta <- prepare_data(meta)



load(file = "data/Lasso80train.data.r")
load(file = "data/Lasso80test.data.r")

row.names(train.data) <- gsub("-", "\\.", row.names(train.data))
print(row.names(train.data))
row.names(train.data) <- gsub("MPA", "MPO", row.names(train.data))
row.names(test.data) <- gsub("-", "\\.", row.names(test.data))
print(row.names(test.data))
row.names(test.data) <- gsub("MPA", "MPO", row.names(test.data))

# Extract the feature names from the dataframe
important_features_sorted <- read.table(file = "../Uwes_wishes/Nature_Code/data/important_features_sorted_lasso_21.csv", sep = ",", header =T)
selected_features <- important_features_sorted$Feature  # Assuming the column 'Feature' contains the feature names
# Add the response variable 'activ_vs_remission' to the selected feature list
selected_features <- c(selected_features, "group")
train.data_subset_p <- train.data[, colnames(train.data) %in% colnames(lasso_global)]
train.data_subset_p <- prepare_data(train.data_subset_p)
train.data_subset_p$Experimental_Group <- NULL
test.data_subset_p <- test.data[, colnames(test.data) %in% colnames(lasso_global)]
test.data_subset_p <- prepare_data(test.data_subset_p)
test.data_subset_p$Experimental_Group <- NULL

# Lasso Features for Evaluation
formula.final <- group ~ LRG1 + B2M.1 + AHSG + CLEC3B + COMP + IGFBP3 + MCAM + TNXB +C9 + F9 + MRC1 + PKM + APMAP + CDH2 + HSPA8 + KRT2 + KRT6B + H6PD + HSPA8 + KRT78
# Convert group to factor with explicit levels
train.data_subset_p$group <- factor(train.data_subset_p$group, levels = c(0, 1))
test.data_subset_p$group <- factor(test.data_subset_p$group, levels = c(0, 1))

# Train the  model
set.seed(123) # For reproducibility
model <- glm(formula.final, data = train.data_subset_p, family = binomial())

# predict on test
test_probabilities <- predict(model, newdata = test.data_subset_p, type = "response")
test_predictions <- ifelse(test_probabilities > 0.5, 1, 0)
test_predictions <- factor(test_predictions, levels = c(0,1))

# Confusion Matrix
conf_matrix <- confusionMatrix(test_predictions, test.data_subset_p$group)

# Performance-Matrix
precision <- conf_matrix$byClass["Pos Pred Value"]
recall <- conf_matrix$byClass["Sensitivity"]
specificity <- conf_matrix$byClass["Specificity"]
npv <- conf_matrix$byClass["Neg Pred Value"]
f1_score <- 2 * (precision * recall) / (precision + recall)


actual_outcomes <- as.numeric(as.character(test.data_subset_p$group))
roc_obj <- roc(actual_outcomes, test_probabilities)

results_df <- data.frame(
  Metric = c("AUC", "Precision", "Recall", "Specificity", "NPV", "F1 Score"),
  Value = c(auc(roc_obj), precision, recall, specificity, npv, f1_score)
)

openxlsx::write.xlsx(list(
  Confusion_Matrix = as.data.frame.matrix(conf_matrix$table),
  Performance_Metrics = results_df,
  ROC_Data = data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  )
), "../Uwes_wishes/data/glmvon21_matrix_evaluation_final_test_results.xlsx")


plot(roc_obj,
     col = "blue1",
     lwd = 2,
     main = "ROC Curve for 21er Panel (glm)",
     print.auc = TRUE,
     print.auc.col = "blue1")
abline(a = 0, b = 1, lty = 2, col = "gray")
legend("bottomright", legend = c("21er Panel (glm)"), col = "blue1", lwd = 2, bty = "n")

# Uwe Plot for Confusionmatrix
all_samples_df <- data.frame(
  Sample = rownames(test.data_subset_p),
  Predicted = test_predictions,
  Actual = test.data_subset_p$group,
  Probability = test_probabilities
)

all_class <- all_samples_df %>%
  mutate(
    Actual = as.numeric(as.character(Actual)),
    Predicted = as.numeric(as.character(Predicted)),
    Jittered_Actual = Actual + runif(n(), -0.4, 0.4),
    Jittered_Predicted = Predicted + runif(n(), -0.4, 0.4)
  )

uwe_plot <- ggplot(all_class, aes(x = Jittered_Actual, y = Jittered_Predicted, fill = as.factor(Predicted))) +
  geom_point(size = 4, alpha = 0.8, pch = 21, color = "black") +
  scale_fill_manual(values = c("brown", "#708238")) +
  scale_x_continuous(breaks = c(0, 1), labels = c("Active", "Remission"), limits = c(-0.5, 1.5)) +
  scale_y_continuous(breaks = c(0, 1), labels = c("Active", "Remission"), limits = c(-0.5, 1.5)) +
  labs(
    x = "Actual Class",
    y = "Predicted Class",
    fill = "Prediction"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    panel.grid = element_blank()
  )

# ggsave(uwe_plot, file = "../Uwes_wishes/April_Figure/21_panelconfusionMatrix_glm_UweStyle.svg", width = 5, height = 4)

# collinearty == correlation
df.rf <- read.table(file = "data/605_data_unique_genenames.tsv", sep ="\t", header = T, row.names = 1)
# Load metadata
meta <- read.table(file = "data/all_metadata.tsv", sep = "\t", header = T, row.names = 1)
# Data preprocessing
meta <- meta %>%
  filter(!group=="HC")
meta$Pseudonym <- NULL
meta$Sample.Number <- NULL
meta <- meta %>%
  mutate(
    sex = if_else(Sex == 'female', 0, 1),
    activ_vs_remission = if_else(disease.nr == "active", 0, 1)
  )
# Additional data preprocessing
row.names(meta) <- meta$Sample_ID_Proteomics
meta$Sample_ID_Proteomics <- NULL
meta$group <- NULL
meta$ANCA_specifity_Disease <- NULL
meta$disease.nr <- NULL

important_features_sorted <- read.table(file = "data/important_features_sorted_lasso_21.csv", sep = ",", header =T)

library(corrplot)
# Extract only the 21 important features from df.rf
important_features <- important_features_sorted$Features # Assuming the column name is 'Feature_Name'
#important_features <- important_features[-1, ]
df.rf_filtered <- df.rf[, colnames(df.rf) %in% important_features]

# Merge data with metadata based on row names
combined_data <- merge(df.rf_filtered, meta, by = "row.names")
row.names(combined_data) <- combined_data$Row.names
combined_data$Row.names <- NULL
combined_data <- combined_data[, c(1:21)]

# Calculate correlations between features and metadata columns
cor_results <- cor(combined_data, use = "pairwise.complete.obs", method = "spearman")

# View correlation matrix
print(cor_results)

# Generate correlation plot
cor_mat_1dec <- round(cor_results, 1)


corrplot(
  cor_results,
  method       = "color",
  type         = "lower",
  tl.col       = "black",
  tl.srt       = 45,
  tl.cex       = 1.2,                                             
  addCoef.col  = NULL,                                            
  cl.cex       = 1.1,                                             
  col          = colorRampPalette(c("lightgrey", "white", "blue"))(400), 
  mar          = c(1, 1, 1, 1)                                    
)

corrplot(
  cor_results,
  method       = "color",
  type         = "lower",
  diag         = FALSE,    
  tl.col       = "black",
  tl.srt       = 45,
  tl.cex       = 1.2,
  addCoef.col  = NULL,
  cl.cex       = 1.1,
  col          = colorRampPalette(c("red", "white", "blue"))(400),
  mar          = c(1,1,1,1)
)


corrplot(
  cor_results,
  method       = "color",                            # Farb‐Tiles
  type         = "lower",                            # Unteres Dreieck
  tl.col       = "black",                            # Label‐Farbe
  tl.srt       = 45,                                 # Label‐Rotation
  tl.cex       = 1.2,                                # Label‐Schriftgröße
  addCoef.col  = "black",                            # Zahlen in Schwarz
  number.cex   = 0.9,                                # Zahlen‐Schriftgröße etwas kleiner
  number.digits= 1,                                  # Eine Nachkommastelle
  cl.cex       = 1.1,                                # Legenden‐Schriftgröße
  cl.align.text= "l",                                # Legendentext linksbündig
  col          = colorRampPalette(c("red","white","blue"))(400),
  mar          = c(1, 1, 1, 1)                       # Kleine Ränder
)



# 0.5 - 0.7: Moderate to strong correlation

# Create a mask for high correlations
high_cor <- abs(cor_results) > 0.8 & abs(cor_results) < 1
# Get the row and column indices of all TRUE values in high_cor
high_cor_indices <- which(high_cor, arr.ind = TRUE)

# high_cor_indices
# row col
# KRT6B  18  17
# KRT78  21  17
# KRT2   17  18
# KRT2   17  21


# PRM semi pure
lasso_global <- read.csv("data/21_lasso_global.csv", row.names=1)
# load the PRM data
prm <- read.table(file ="data/rpm_wo_contols.csv", sep =",", header =T, row.names = 1) 

prm$Experimental_Group <- row.names(prm)
prm <-  prm %>%
  mutate(group = case_when(
    grepl("active", Experimental_Group) ~ "active",
    grepl("remission", Experimental_Group) ~ "remission",
    grepl("remisson", Experimental_Group) ~ "remission",
    grepl("Remission", Experimental_Group) ~ "remission",
    grepl("HC", Experimental_Group) ~ "healthy",
    TRUE ~ NA_character_ # for words that do not match any of the conditions
  ))

prm <- prm %>% 
  filter(!group == "healthy")


df.red <- prm%>% 
  mutate(group = if_else(group == "active", 0,1)) # active = 0
## active == 1, remission == 0
df.red$Experimental_Group <- NULL
df.red <- df.red %>% 
  drop_na()



# Aggregates df.red (PRM peptide data) to protein level
# Standardizes sample names
# Merges with lasso_global based on sample names
# Calculates protein-wise correlations between PRM and LASSO expression
# outputs a correlation summary and an optional heatmap

df.red$Sample <- rownames(df.red)

# Convert from wide to long format
df_long <- df.red %>%
  pivot_longer(cols = -Sample, names_to = "Peptide", values_to = "Value") %>%
  mutate(Protein = sub("_.*", "", Peptide))  # Extract protein name before underscore

# Aggregate peptide intensities to protein level (mean per sample × protein)
df_protein <- df_long %>%
  group_by(Sample, Protein) %>%
  summarise(Mean_Intensity = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Protein, values_from = Mean_Intensity)

df_protein <- df_protein %>%
  mutate(Sample_clean = Sample) %>%
  
  # Step 1: PRM name cleanup
  mutate(
    Sample_clean = gsub("\\.", "_", Sample_clean),               # first: replace . with _
    Sample_clean = gsub("-", "_", Sample_clean),                 # then: ensure separator is _
    
    # Step 2: Fix disease type mapping
    Sample_clean = gsub("PR3_active", "PR3_active", Sample_clean),
    Sample_clean = gsub("PR3_remission", "PR3_Remission", Sample_clean),
    Sample_clean = gsub("PR3_Remission", "PR3_Remission", Sample_clean),
    Sample_clean = gsub("MPA_active", "MPO_active", Sample_clean),
    Sample_clean = gsub("MPA_Remission", "MPO_Remission", Sample_clean),
    
    # Step 3: Move disease label to front: A20_26_PR3_active → PR3_active_A20.26
    Sample_clean = sub("^(.*)_(PR3|MPO)_(active|Remission)$", "\\2_\\3_\\1", Sample_clean),
    
    # Step 4: Convert underscores in ID back to dots (e.g. A20_26 → A20.26)
    Sample_clean = gsub("([A-Z][0-9]+)_([0-9]+)", "\\1.\\2", Sample_clean)
  )
lasso_global$Sample <- rownames(lasso_global)
setdiff(lasso_global$Sample, df_protein$Sample_clean)
setdiff(df_protein$Sample_clean, lasso_global$Sample)

# Merge on cleaned sample names
merged_data <- inner_join(df_protein, lasso_global, by = c("Sample_clean" = "Sample"))

# Get the matched PRM and LASSO columns
prm_cols <- grep("\\.x$", names(merged_data), value = TRUE)
lasso_cols <- grep("\\.y$", names(merged_data), value = TRUE)

# Optional: align order by protein name
clean_names <- function(x) gsub("\\.x$|\\.y$", "", x)
common_proteins <- intersect(clean_names(prm_cols), clean_names(lasso_cols))

# Reorder both sets
prm_cols <- paste0(common_proteins, ".x")
lasso_cols <- paste0(common_proteins, ".y")

# Calculate correlation per patient
per_patient_cor <- data.frame(
  Sample = merged_data$Sample_clean,
  Pearson_r = mapply(function(i) {
    cor(
      as.numeric(merged_data[i, prm_cols]),
      as.numeric(merged_data[i, lasso_cols]),
      use = "pairwise.complete.obs"
    )
  }, 1:nrow(merged_data))
)

# View results
print(per_patient_cor)

# Build matrices
prm_mat <- as.matrix(merged_data[, prm_cols])
lasso_mat <- as.matrix(merged_data[, lasso_cols])

# Rename columns for heatmap
colnames(prm_mat) <- common_proteins
colnames(lasso_mat) <- common_proteins

# clean names
rownames(protein_cormat) <- gsub("PRM: ", "", rownames(protein_cormat))
colnames(protein_cormat) <- gsub("Proteome: ", "", colnames(protein_cormat))

ordered_names <- sort(intersect(rownames(protein_cormat), colnames(protein_cormat)))
protein_cormat_ordered <- protein_cormat[ordered_names, ordered_names]

# Add column/row labelling again
rownames(protein_cormat_ordered) <- paste0("PRM: ", ordered_names)
colnames(protein_cormat_ordered) <- paste0("Proteome: ", ordered_names)

# Heatmap 
pheatmap(
  protein_cormat_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
)

prm_proteins <- unique(sub("_.*", "", colnames(df.red)))

# Initialise results list
all_correlation_results <- list()

# itterate for each Protein
for (target_protein in prm_proteins) {
  
  peptides <- grep(paste0("^", target_protein, "_"), colnames(df.red), value = TRUE)
  
  if (length(peptides) == 0 || !(target_protein %in% colnames(lasso_global))) {
    next  # skip if protein or LASSO- protein is not available oder
  }
  
  # PRM-values + Sample-Namen
  prm_values <- df.red[, peptides, drop = FALSE]
  prm_values$Sample <- rownames(df.red)
  
  # LASSO-values
  global_values <- lasso_global[, target_protein]
  names(global_values) <- rownames(lasso_global)
  
  prm_values <- prm_values %>%
    mutate(Sample_clean = Sample) %>%
    mutate(
      Sample_clean = gsub("^X", "", Sample_clean),
      Sample_clean = gsub("\\.", "_", Sample_clean),
      Sample_clean = gsub("-", "_", Sample_clean),
      Sample_clean = gsub("PR3_active", "PR3_active", Sample_clean),
      Sample_clean = gsub("PR3_remission", "PR3_Remission", Sample_clean),
      Sample_clean = gsub("MPA_active", "MPO_active", Sample_clean),
      Sample_clean = gsub("MPA_Remission", "MPO_Remission", Sample_clean),
      Sample_clean = sub("^(.*)_(PR3|MPO)_(active|Remission)$", "\\2_\\3_\\1", Sample_clean),
      Sample_clean = gsub("([A-Z][0-9]+)_([0-9]+)", "\\1.\\2", Sample_clean)
    )
  
  # Matchen
  matched <- intersect(prm_values$Sample_clean, rownames(lasso_global))
  
  if (length(matched) < 3) next on
  
  # calculate Pearson r
  cor_peptides <- sapply(peptides, function(pept) {
    x <- prm_values[prm_values$Sample_clean %in% matched, pept]
    y <- global_values[matched]
    
    valid <- complete.cases(x, y)
    if (sum(valid) >= 2) {
      cor(x[valid], y[valid], use = "complete.obs")
    } else {
      NA
    }
  })
  
  best_idx <- which.max(abs(cor_peptides))
  cor_df <- data.frame(
    Protein = target_protein,
    Peptide = names(cor_peptides),
    Pearson_r = round(cor_peptides, 3)
  )
  
  all_correlation_results[[target_protein]] <- cor_df
}

# Combine everything into one DataFrame
all_peptides_correlation_df <- do.call(rbind, all_correlation_results)

#all peptides per Protein
print(all_peptides_correlation_df)

ggplot(all_peptides_correlation_df, aes(x = reorder(Peptide, Pearson_r), y = Pearson_r, fill = Protein)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Correlation of PRM Peptides with Global Proteomics (LASSO)",
    x = "Peptide",
    y = "Pearson Correlation (r)"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 12),
        axis.titel.x = element_blank()) -> prearsonCorr

# ggsave(prearsonCorr, file = "../Uwes_wishes/April_Figure/prearsonCorr_prm_proteomics.svg", width = 9, height = 12)


# In wide format 
heatmap_df <- all_peptides_correlation_df %>%
  pivot_wider(names_from = Peptide, values_from = Pearson_r)

# rownames replace
heatmap_mat <- as.data.frame(heatmap_df)
rownames(heatmap_mat) <- heatmap_mat$Protein
heatmap_mat$Protein <- NULL
heatmap_mat <- as.matrix(heatmap_mat)

pheatmap(
  heatmap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  na_col = "grey90",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Peptide-Protein Correlation (Pearson r)"
)


heatmap_long <- as.data.frame(as.table(heatmap_mat)) %>%
  rename(Protein = Var1, Peptide = Var2, Correlation = Freq)

# remove NA
heatmap_long <- heatmap_long %>% filter(!is.na(Correlation))

# Optional: Show only lower half (if protein == peptide is possible, otherwise omit)
# Here: all against all → no triangular matrix, but you can specify a ‘step-like order’
heatmap_long <- heatmap_long %>%
  mutate(Protein = factor(Protein, levels = unique(Protein)),
         Peptide = factor(Peptide, levels = unique(Peptide)))

# Plot
ggplot(heatmap_long, aes(x = Peptide, y = Protein, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limits = c(-1, 1), na.value = "grey90") +
  labs(
    title = "Correlation between PRM peptides and global Proteins",
    x = "PRM Peptides",
    y = "Proteins",
    fill = "Pearson r"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) -> heat

# ggsave(heat, file = "../Uwes_wishes/April_Figure/prearsonCorr_prm_proteomicsheat.svg", width = 15, height = 12)


the15 <-read.csv2("../Uwes_wishes/Nature_Code/data/Panel_21_PRM_20250115_FC_pvalue (3).csv",na.strings = "")


row.names(the15) <- the15$Peptide
the15$Peptide_label <- gsub("_[0-9]+$", "", rownames(the15))

# Assign category correctly
the15$Category <- "Other"
the15$Category[the15$negative.control == " Yes"] <- "Negative Control"
the15$Category[the15$postive.control == " Yes"] <- "Positive Control"
the15$Category[the15$panel.21 == "Yes"] <- "Panel 21"
the15$Category[the15$panel.7.peptides == "Yes"] <- "Panel 7 Peptide"

# binär significance
the15$Significance <- ifelse(
  !is.na(the15$Welch.s.T.test.Significant.active_remission_FDR5) &
    the15$Welch.s.T.test.Significant.active_remission_FDR5 == "Yes",
  "sig",
  "nonsig"
)
# groups per color
the15$fill_group <- paste0(the15$Category, "_", the15$Significance)

# Shape: only panel 7 peptides filled (16), rest empty (21)
the15$Shape <- ifelse(the15$Category == "Panel 7 Peptide", 16, 21)


# Get Panel 7 and non-Panel 7
filled_points <- subset(the15, Category == "Panel 7 Peptide")
outline_points <- subset(the15, Category != "Panel 7 Peptide")

# Define colours for outline circles (coloured border, filled white)
outline_colors <- c(
  "Panel 21_sig" = "#238b45",
  "Panel 21_nonsig" = "#a1d99b",
  "Positive Control_sig" = "#08519c",
  "Positive Control_nonsig" = "#9ecae1",
  "Negative Control_sig" = "gray40",
  "Negative Control_nonsig" = "gray80",
  "Other_sig" = "black",
  "Other_nonsig" = "gray90"
)

# Colours for panel 7 filled dots
filled_colors <- c(
  "Panel 7 Peptide_sig" = "#c70000",# dark red
  "Panel 7 Peptide_nonsig" = "#ffada3" # light red
)

# Total colours for scale_fill_manual (for legend)
fill_legend <- c(outline_colors, filled_colors)


ggplot() +
  # Outline points (empty circles with colored boader)
  geom_point(
    data = outline_points,
    aes(
      x = logFC.active_remission,
      y = p.value.active_remission...log10.,
      color = fill_group
    ),
    fill = "white", shape = 21, size = 4, stroke = 1, alpha = 0.9
  ) +
  # Panel 7: colored points
  geom_point(
    data = filled_points,
    aes(
      x = logFC.active_remission,
      y = p.value.active_remission...log10.,
      fill = fill_group
    ),
    shape = 21, color = "black", size = 4, stroke = 0.5, alpha = 0.9
  ) +
  # Labels for Panel 7
  geom_label_repel(
    data = filled_points,
    aes(
      x = logFC.active_remission,
      y = p.value.active_remission...log10.,
      label = Peptide
    ),
    size = 3.5,
    fill = "white",
    box.padding = 0.25,
    label.padding = 0.15,
    label.size = 0.3,
    segment.color = "gray30",
    max.overlaps = 20
  ) +
  scale_color_manual(
    values = outline_colors,
    guide = "legend",
    name = "Peptide Category"
  ) +
  scale_fill_manual(
    values = filled_colors,
    guide = "legend",
    name = "Peptide Category"
  ) +
  labs(
    x = "Fold Change (Active vs. Remission)",
    y = "-log10(p-value)"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )


volcano_plot <- ggplot(the15, aes(
  x = logFC.active_remission,
  y = p.value.active_remission...log10.,
fill = Category,
  color = Welch.s.T.test.Significant.active_remission_FDR5
)) +
  geom_point(
    shape = 21, size = 4, stroke = 0.8, alpha = 0.8
  ) +
  geom_label_repel(
    data = subset(the15, Category == "Panel 7 Peptide"),
    aes(label = Peptide),
    size = 3.5,
    fill = "white",       
    box.padding = 0.25,
    label.padding = 0.15,
    label.size = 0.3,
    segment.color = "gray30",
    max.overlaps = 20
  ) +
  scale_color_manual(
    values = c(
      "Negative Control" = "#F2f0e6",
      "Positive Control" = "#BDe7bd",
      "Panel 7 Peptide" = "black",
      "Other" = "#F2f0e6"
    )
    )+
  scale_fill_manual(
    values = c(
      "Negative Control" = "white",
      "Positive Control" = "#BDe7bd",
      "Panel 7 Peptide" = "#7d002b",
      "Other" = "white"
    )
  ) +
  scale_color_manual(
    name = "Significant (FDR < 0.05)",
    values = c("Yes" = "black", "No" = "gray70", "NA" = "gray70"),
    na.value = "gray70"
  ) +
  labs(
    x = "Log2 Fold Change (Active vs. Remission)",
    y = "-log10(p-value)",
    fill = "Peptide Category"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

volcano_plot

the15 <-read.csv2("data/Panel_21_PRM_20250115_FC_pvalue (3).csv",na.strings = "")

the15$Peptide <- gsub(" ", "", the15$Peptide)
# remove empty peptide names (‘’ after trim)
the15 <- the15[the15$Peptide != "", ]
# Optional: adjusted label version for labelling
the15$Peptide_label <- gsub("_[0-9]+.*$", "", the15$Peptide)

the15$Category <- "Other"
the15$Category[the15$negative.control == " Yes"] <- "Negative Control"
the15$Category[the15$postive.control == " Yes"] <- "Positive Control"
the15$Category[the15$panel.7.peptides == "Yes"] <- "Panel 7 Peptide"
the15$Peptide <- gsub("_", " ",the15$Peptide )

volcano_plot <- ggplot(the15, aes(
  x = logFC.active_remission,
  y = p.value.active_remission...log10.
)) +
  # All other peptides (not Panel 7 and not Positive Control)
  geom_point(
    data = subset(the15, !(Category %in% c("Panel 7 Peptide", "Positive Control"))),
    aes(color = Welch.s.T.test.Significant.active_remission_FDR5),
    shape = 21,
    size = 4,
    stroke = 0.8,
    fill = NA,
    alpha = 0.8
  ) +
  # Positive controls: green border only, no fill
  geom_point(
    data = subset(the15, Category == "Positive Control"),
    aes(x = logFC.active_remission, y = p.value.active_remission...log10.),
    shape = 21,
    size = 4,
    stroke = 0.8,
    color = "#00B050",  # Grün
    fill = NA,
    alpha = 0.8
  ) +
  # Panel 7 Peptides: filled in red, border according to significance
  geom_point(
    data = subset(the15, Category == "Panel 7 Peptide"),
    aes(color = Welch.s.T.test.Significant.active_remission_FDR5),
    shape = 21,
    size = 4,
    stroke = 0.8,
    fill = "#c70000",
    alpha = 0.8
  ) +
  # Label panel 7 peptides
  geom_label_repel(
    data = subset(the15, Category == "Panel 7 Peptide"),
    aes(label = Peptide_label),
    size = 3.5,
    fill = "white",
    box.padding = 0.25,
    label.padding = 0.15,
    label.size = 0.3,
    segment.color = "gray30",
    max.overlaps = 20
  ) +
  # Colour legend for significance
  scale_color_manual(
    name = "Significant (FDR < 0.05)",
    values = c("Yes" = "black", "No" = "gray70", "NA" = "gray70"),
    na.value = "gray70"
  ) +
  labs(
    x = "Log2 Fold Change AAV-A vs. AAV-R",
    y = "-log10(p-value)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14)
  )

print(volcano_plot)



boruta_df <- read.csv2(file = "Uwes_wishes/Nature_Code/data/boruta_importance_decision_dataframe.csv")
# Add a column for color based on boruta decision
boruta_df$Color <- ifelse(boruta_df$Decision == "Confirmed", "#c70000",
                          ifelse(boruta_df$Decision == "Tentative",  "lightgrey", "#595959"))

# sort dataframe for importance
boruta_df <- boruta_df[order(boruta_df$Importance, decreasing = TRUE), ]
boruta_df$Feature <- gsub("_", " ",boruta_df$Feature )
boruta_df$Feature <- gsub("\\.", "",boruta_df$Feature )
boruta_df <- boruta_df[boruta_df$Feature != "", ]
boruta_df$Feature <- gsub(" [0-9]+.*$", "", boruta_df$Feature)

ggplot(boruta_df, aes(x = reorder(Feature, Importance), y = Importance, fill = Color, colour = "black")) +
  geom_bar(stat = "identity", colour = "black") +
  coord_flip() + # Für horizontale Balken
  scale_fill_manual(values = c("#c70000" = "#c70000", "lightgrey" = "lightgrey",
                               "#595959" = "#595959"),
                    name = "Boruta Decision") +
  labs(
       x = "Feature",
       y = "Mean Importance (Random Forest Accuracy Decrease)") +
  theme_bw()+
  theme(
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank()
  )

# Boruta selected features
selected_features <- c("MCAM_047", "AHSG_004", "LRG1_042", "CLEC3B_020", "COMP_023", "MRC1_073", "F9_027")

# Add column to mark the dataset
train_data_long <- PRM_semi_pure_lassofeatures80iger_split %>%
  select(all_of(selected_features), group) %>%
  pivot_longer(cols = -group, names_to = "Protein", values_to = "Value") %>%
  mutate(Dataset = "Train")

test_data_long <- PRM_semi_pure_lassofeatures20iger_split %>%
  select(all_of(selected_features), group) %>%
  pivot_longer(cols = -group, names_to = "Protein", values_to = "Value") %>%
  mutate(Dataset = "Test")

# combine training and test data
combined_data_long <- bind_rows(train_data_long, test_data_long)

# Transform 'group'-variable in factor
combined_data_long$group <- factor(combined_data_long$group)

ggplot(combined_data_long, aes(x = Protein, y = Value, color = group)) +
  geom_boxplot(aes(fill = group), alpha = 0.7, outlier.shape = NA, width = 0.8) + # Boxplots etwas breiter
  geom_point(data = subset(combined_data_long, Dataset == "Train"), color = "grey",
             position = position_jitter(width = 0.2), size = 1, alpha = 0.5) + # Transparenz für Trainingsdatenpunkte
  geom_point(data = subset(combined_data_long, Dataset == "Test"), color = "blue",
             position = position_jitter(width = 0.2), size = 1, alpha = 0.5) + # Transparenz für Testdatenpunkte
  scale_fill_manual(values = c("1" = "darkgreen", "0" = "brown"), name = "Klinischer Status") +
  scale_color_manual(values = c("1" = "darkgreen", "0" = "brown"), guide = "none") +
  labs(x = "Protein",
       y = "log2(light/heavy)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Set the order of the levels for the dataset
combined_data_long$Dataset <- factor(combined_data_long$Dataset, levels = c("Test", "Train"))
combined_data_long$Protein <- gsub("_", " ", combined_data_long$Protein)
combined_data_long$Protein <- gsub("\\.", "", combined_data_long$Protein)


ggplot(combined_data_long, aes(x = Protein, y = Value, fill = group,
                               group = interaction(Protein, Dataset, group))) +
  geom_boxplot(aes(color = Dataset),
               position = position_dodge(width = 2),
               width = 2,
               outlier.shape = NA,
               alpha = 0.7) +
  scale_color_manual(values = c("Test" = "blue", "Train" = "grey"), guide = "none") +
  scale_fill_manual(values = c("1" = "darkgreen", "0" = "brown"),
                    name = "clinical status") +
  labs(x = "Protein",
       y = "log2(light/heavy)") +
  theme_bw()+
  theme(
    legend.position = "right",
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank()
  )

# Calculate the median for each combination of protein, dataset and group in order to place the text sensibly
summary_data <- combined_data_long %>%
  group_by(Protein, Dataset, group) %>%
  summarize(med = median(Value), .groups = "drop")

ggplot(combined_data_long, aes(x = Protein, y = Value, fill = group,
                               group = interaction(Protein, Dataset, group))) +
  geom_boxplot(aes(color = Dataset),
               position = position_dodge(width = 1.1),
               width = 1,
               outlier.shape = NA,
               alpha = 0.7) +
  # Textannotation: Zeigt "Train" oder "Test" oberhalb des jeweiligen Boxplots an
  scale_color_manual(values = c("Test" = "blue", "Train" = "grey"), guide = "none") +
  scale_fill_manual(values = c("1" = "darkgreen", "0" = "brown"),
                    name = "clinical status") +
  labs(x = "Protein",
       y = "log2(light/heavy)") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank())

combined_data_long$Protein <- gsub("_", " ", combined_data_long$Protein)
combined_data_long$Protein <- gsub("\\.", "", combined_data_long$Protein)


ggplot(combined_data_long, aes(x = Value, y = Protein, fill = Dataset)) +
  # violin plot shows the density of the distribution per protein.
  geom_violin(alpha = 0.6, trim = FALSE, 
              scale = "width", 
              position = position_dodge(width = 0.9)) +
  # narrow boxplot showing central values (median, quartiles).
  geom_boxplot(width = 0.2, outlier.shape = NA, 
               alpha = 0.8, 
               position = position_dodge(width = 0.9)) +
  # mean value as an additional marker
  stat_summary(fun = mean, geom = "point", shape = 21, fill = "white",
               size = 3, position = position_dodge(width = 0.9)) +
  # faceting by clinical status (group) - separate panels for e.g. 0 and 1.
  facet_wrap(~ group, scales = "free_y") +
  scale_fill_manual(values = c("Train" = "grey", "Test" = "blue"), 
                    name = "Dataset") +
  labs(x = "log2(light/heavy)")+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.position = "right")



combined_data_long <- combined_data_long %>%
  mutate(group_label = ifelse(group == "0", "active AAV", "rem AAV"))

ggplot(combined_data_long, aes(x = Value, y = Protein)) +
# Move violin ‘Train’ (pseudo half) slightly upwards 
geom_violin(
  data = subset(combined_data_long, Dataset == "Train"),
  aes(fill = Dataset),
  alpha    = 0.6,
  scale    = "width",
  trim     = FALSE,
  position = position_nudge(y =  0.2)) +
# Move violin ‘Test’ (pseudo half) slightly downwards
geom_violin(
  data = subset(combined_data_long, Dataset == "Test"),
  aes(fill = Dataset),
  alpha    = 0.6,
  scale    = "width",
  trim     = FALSE,
  position = position_nudge(y = -0.2)
) +
# Show all points in the ‘centre’
geom_point(
  aes(color = Dataset),
  position = position_jitter(width = 0, height = 0),
  alpha    = 0.7,
  size     = 2
) +
  facet_wrap(~ group_label, scales = "free_y") +
  scale_fill_manual(values = c("Train" = "grey", "Test" = "blue")) +
  labs(x = "log2(light/heavy)", y = "Protein") +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 14),
    strip.text = element_text(size = 16)
  )


# We need the following columns present in the data frame ‘combined_data_long’:
# - Protein: Name/ID of the protein
# - Value: e.g. log2(light/heavy) values
# - Dataset: with the values ‘Train’ (80%) and ‘Test’ (20%)
# - group_label: e.g. ‘active’ and ‘remission’

# Create new combination variable (for separate X-axis categories)
combined_data_long <- combined_data_long %>%
  mutate(ProtGroup = interaction(Protein, group_label, sep = " - "))
combined_data_long$Protein <- gsub(" [0-9]+.*$", "", combined_data_long$Protein)

# wilcox-test
wilcox_results <- combined_data_long %>%
  group_by(Protein) %>%
  wilcox_test(Value ~ group_label) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  mutate(
    p.adj.label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ paste0("p=", round(p.adj, 3))
    )
  )

ggplot(
  combined_data_long, 
  aes(
    x    = Protein,
    y    = Value,
    fill = group_label
  )
) +
  # Train-Data (left half)
  geom_half_violin(
    data     = subset(combined_data_long, Dataset == "Train"),
    side     = "l",draw_quantiles = 0.5,
    width    = 2,          
    scale    = "width",      
    trim     = FALSE,       
    alpha    = 0.5,
    position = position_dodge(width = 1.2)  
  ) +
  # Test-Data (right half)
  geom_half_violin(
    data     = subset(combined_data_long, Dataset == "Test"),
    side     = "r",draw_quantiles = 0.5,
    width    = 2,
    scale    = "width",
    trim     = FALSE,
    alpha    = 0.5,
    position = position_dodge(width = 1.2)
  ) +
  # Punkte mit Jitter
  geom_point(
    aes(fill = group_label),
    pch = 21,
    position = position_jitterdodge(
      dodge.width  = 1.2,
      jitter.width = 0.3
    ),
    size  = 3,
    alpha = 0.8,
    color = "black"
  )+
  # Signifikanz-Annotationen
  geom_text(
    data = wilcox_results,
    aes(
      x = Protein,
      y = max(combined_data_long$Value) * 1.15,
      label = p.adj.label
    ),
    inherit.aes = FALSE,
    size = 5,
    vjust = -0.5
  ) +
  # Farben und Theme
  scale_fill_manual(
    values =c("active AAV" = "#ffd973", "rem AAV" = "#0000ab"),
    name   = "disease status"
  ) +
  labs(
    x = NULL,
    y = "log2(light/heavy)",
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(
      size = 12,
    ),
    axis.text.y = element_text(
      size = 12
    ),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(min(combined_data_long$Value), 
                           max(combined_data_long$Value) * 1.1))+
  coord_flip()



Marie_log2_welchT7 <- read.csv2("../Uwes_wishes/data/Panel_7_PRM_TQL.csv")
# Create the volcano plot

# panel proteins
proteins <- c("AHSG_004","CLEC3B_020", "COMP_023", "F9_025", "LRG1_045", "MCAM_048", "MRC1_073")

# Assuming 'names_vector' contains your target inflammatory genes and is correctly formatted
Marie_log2_welchT7 <- Marie_log2_welchT7 %>%
  mutate(
    protein = ifelse(Peptide %in% proteins, "Panel", "nonPanel")
  )



volcano_plot <- ggplot(Marie_log2_welchT7, aes(
  x = Difference.Active_Remission..log2.FC.,
  y = p.value.Active_Remission...log10.
)) +
  geom_point(
    data = subset(Marie_log2_welchT7, protein == "Panel"),
    aes(color = Welch.s.T.test.Significant.Prag_Active_Prag_Remission),
    size = 3,  # Bigger size for Inflammatory
    shape = 21,  # Using circles for both, can change if needed
    fill = "darkred",  # Fill red for inflammatory
    stroke = 1, alpha = 0.8
  ) +
  geom_point(
    data = subset(Marie_log2_welchT7, protein == "nonPanel"),
    aes(color = Welch.s.T.test.Significant.Prag_Active_Prag_Remission),
    size = 2,  # Slightly smaller size for Non-inflammatory
    shape = 21,  # Using circles for both, can change if needed
    fill = "darkgreen",  # Fill dark yellow for non-inflammatory
    stroke = 1, alpha = 0.5
  ) +
  geom_text_repel(
    aes(label = ifelse(protein == "Panel", Peptide, "")),
    max.overlaps = 10, size = 3
  )+
  scale_color_manual(values = c("significant" = "black", "non_significant" = "gray")) +
  labs(
    x = "Log2 Fold Change active vs remission",
    y = "-log10(p-value)",
    color = "Significance",
    fill = "Protein Type"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

print(volcano_plot)

ggsave(volcano_plot, file = "../Uwes_wishes/April_Figure/2549_volcanos7ner.svg", width = 5, height = 6)