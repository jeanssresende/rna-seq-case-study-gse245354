# Lesson 08 – PCA and Heatmap Analysis

# ==========================================================
#
# RNA-Seq Case Study
# Dataset: GSE245354
#
# Lesson 08
# PCA and Heatmap Analysis
#
# Script:
# 08_pca_and_heatmap_analysis.R
#
# Author: Jean Resende
#
# ==========================================================


# ==========================================================
# Lesson Objectives
# ==========================================================

# In this lesson you will learn:
#
# • Why multivariate analysis is important
# • How to perform Principal Component Analysis (PCA)
# • How to visualize sample relationships
# • How to identify highly variable genes
# • How to build a heatmap
# • How to interpret clustering patterns
#
# ==========================================================


# ==========================================================
# Installing required packages
# ==========================================================

if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")


# ==========================================================
# Loading packages
# ==========================================================

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(stringr)


# ==========================================================
# Loading expression matrix
# ==========================================================

counts <- read.delim(
  "data/raw/GSE245354/GSE245354_Meegan_Normalized_counts.txt.gz"
)


# ==========================================================
# Loading metadata
# ==========================================================

metadata <- read.csv(
  "data/metadata/sample_metadata.csv"
)


# ==========================================================
# Preparing expression matrix
# ==========================================================

expr <- counts[, 9:ncol(counts)]

rownames(expr) <- counts$Gene.Symbol


# Remove genes without gene symbol
expr <- expr[rownames(expr) != "" & !is.na(rownames(expr)), ]


# ==========================================================
# Log2 transformation
# ==========================================================

# RNA-seq data are often highly skewed.
# Log transformation makes distributions more symmetric.

expr_log2 <- log2(expr + 1)


# ==========================================================
# PCA requires samples in rows
# ==========================================================

expr_pca <- t(expr_log2)


# ==========================================================
# Running PCA
# ==========================================================

pca <- prcomp(expr_pca, scale. = TRUE)


# ==========================================================
# Percentage of explained variance
# ==========================================================

variance <- pca$sdev^2

variance_percent <- round(
  variance / sum(variance) * 100,
  2
)

variance_percent[1:5]


# ==========================================================
# Building PCA data frame
# ==========================================================

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)


# ==========================================================
# Preparing metadata for join
# ==========================================================

metadata <- read.csv(
  "data/metadata/sample_metadata.csv",
  stringsAsFactors = FALSE
)


metadata <- metadata %>%
  mutate(
    SampleID = str_extract(Description, "(?<=NP-)\\d+"),
    MatrixSample = paste0("X5958.NP.", SampleID, "_S1")
  )

metadata %>%
  select(SampleID, MatrixSample, Treatment, Time, Donor) %>%
  head()


# ==========================================================
# Integrating PCA coordinates with metadata
# ==========================================================

pca_df <- pca_df %>%
  left_join(
    metadata %>%
      select(MatrixSample, Treatment, Time, Donor),
    by = c("Sample" = "MatrixSample")
  )


# ==========================================================
# Inspecting PCA table
# ==========================================================

head(pca_df)


# ==========================================================
# PCA plot colored by treatment
# ==========================================================

pca_plot_treatment <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Treatment)
) +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA - Colored by Treatment",
    x = paste0("PC1 (", variance_percent[1], "%)"),
    y = paste0("PC2 (", variance_percent[2], "%)")
  )


pca_plot_treatment


# ==========================================================
# PCA plot colored by time
# ==========================================================

pca_plot_time <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Time)
) +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA - Colored by Collection Time",
    x = paste0("PC1 (", variance_percent[1], "%)"),
    y = paste0("PC2 (", variance_percent[2], "%)")
  )


pca_plot_time


# ==========================================================
# PCA plot colored by donor
# ==========================================================

pca_plot_donor <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = factor(Donor))
) +
  geom_point(size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA - Colored by Donor",
    x = paste0("PC1 (", variance_percent[1], "%)"),
    y = paste0("PC2 (", variance_percent[2], "%)"),
    color = "Donor"
  )

pca_plot_donor


# ==========================================================
# Saving PCA plots
# ==========================================================

ggsave(
  "results/pca_treatment.png",
  pca_plot_treatment,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "results/pca_time.png",
  pca_plot_time,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  "results/pca_donor.png",
  pca_plot_donor,
  width = 7,
  height = 5,
  dpi = 300
)


# ==========================================================
# Identifying highly variable genes
# ==========================================================

gene_variance <- apply(expr_log2, 1, var)

gene_variance <- sort(gene_variance, decreasing = TRUE)

head(gene_variance)


# ==========================================================
# Selecting top variable genes
# ==========================================================

top_genes <- names(gene_variance)[1:50]

expr_top <- expr_log2[top_genes, ]


# ==========================================================
# Preparing heatmap annotation
# ==========================================================

annotation_col <- metadata %>%
  select(MatrixSample, Treatment, Time, Donor) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$MatrixSample

annotation_col$MatrixSample <- NULL

# Convert Donor to factor for annotation colors
annotation_col$Donor <- factor(annotation_col$Donor)

# Reorder annotation to match matrix columns
annotation_col <- annotation_col[colnames(expr_top), ]


# ==========================================================
# Heatmap
# ==========================================================

pheatmap(
  expr_top,
  scale = "row",
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(
    rev(brewer.pal(9, "RdBu"))
  )(100),
  show_rownames = FALSE,
  fontsize_col = 8,
  main = "Top 50 Most Variable Genes",
  filename = "results/heatmap_top50_variable_genes.png",
  width = 8,
  height = 10
)


# ==========================================================
# Lesson Summary
# ==========================================================

# Today you learned:
#
# ✓ Log2 transformation
# ✓ Principal Component Analysis
# ✓ Explained variance
# ✓ PCA visualization
# ✓ Gene variance estimation
# ✓ Heatmap construction
# ✓ Hierarchical clustering
# ✓ Sample clustering interpretation
#
# ==========================================================


# ==========================================================
# Thinking Like a Bioinformatician
# ==========================================================

# PCA asks:
#
# "Which samples are globally similar?"
#
# Heatmaps ask:
#
# "Which genes drive those similarities?"
#
# PCA summarizes thousands of genes into a few
# dimensions, while heatmaps display the detailed
# expression patterns of selected genes.
#
# These analyses are exploratory and help detect:
#
# • treatment effects
# • time effects
# • batch effects
# • outlier samples
# • donor-specific patterns
#
# ==========================================================


# ==========================================================
# Challenge
# ==========================================================

# 1. Increase the number of genes in the heatmap
#    from 50 to 100.
#
# 2. Color the PCA by Donor.
#
# 3. Add sample labels to the PCA plot.
#
# 4. Identify samples that cluster together.
#
# 5. Identify genes that show opposite patterns
#    between treatment groups.
#
# ==========================================================
