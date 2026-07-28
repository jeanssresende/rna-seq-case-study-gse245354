# ==========================================================
#
# RNA-Seq Case Study
# Dataset: GSE245354
#
# Lesson 06
# Exploratory Data Analysis (EDA)
#
# Script:
# 06_exploratory_data_analysis.R
#
# Author: Jean Resende
#
# ==========================================================


# ==========================================================
# Lesson Objectives
# ==========================================================

# In this lesson you will learn:
#
# • How to explore an RNA-seq expression matrix
# • How to summarize gene expression values
# • How to visualize expression distributions
# • How to identify highly expressed genes
# • How to identify lowly expressed genes
# • How to inspect sample variability
#
# ==========================================================


# ==========================================================
# Loading packages
# ==========================================================

library(dplyr)


# ==========================================================
# Loading the expression matrix
# ==========================================================

counts <- read.delim("data/raw/GSE245354/GSE245354_Meegan_Normalized_counts.txt.gz")


# ==========================================================
# First look at the data
# ==========================================================
dim(counts)

head(counts)

tail(counts)

names(counts)


# ==========================================================
# Separating annotation from expression values
# ==========================================================

# The first columns contain gene annotations.
# Expression values start after the annotation columns.

expr <- counts[, c(9:ncol(counts))]

dim(expr)


# ==========================================================
# General summary
# ==========================================================

summary(expr)


# ==========================================================
# Minimum and maximum values
# ==========================================================

range(as.matrix(expr))


# ==========================================================
# Expression distribution of one sample
# ==========================================================

hist(
  expr[,1],
  breaks = 50,
  main = "Expression Distribution",
  xlab = "CPM",
  col = "steelblue"
)


# ==========================================================
# Boxplot of all samples
# ==========================================================

boxplot(
  expr,
  las = 2,
  cex.axis = 0.7,
  col = "lightgray",
  main = "Expression Distribution Across Samples",
  ylab = "CPM"
)


# ==========================================================
# Density curve
# ==========================================================

plot(
  density(expr[,1]),
  lwd = 2,
  main = "Density Distribution",
  xlab = "CPM"
)


# ==========================================================
# Mean expression of each gene
# ==========================================================

counts$MeanExpression <- rowMeans(expr)

head(counts)


# ==========================================================
# Highest expressed genes
# ==========================================================

highest_expression <- counts %>%
  arrange(desc(MeanExpression))


head(highest_expression)


# ==========================================================
# Lowest expressed genes
# ==========================================================

lowest_expression <- counts %>%
  arrange(MeanExpression)


head(lowest_expression)


# ==========================================================
# Average expression of each sample
# ==========================================================

sample_means <- colMeans(expr)

sample_means


# ==========================================================
# Sample medians
# ==========================================================

sample_medians <- apply(expr, 2, median)

sample_medians


# ==========================================================
# Sample standard deviation
# ==========================================================

sample_sd <- apply(expr, 2, sd)

sample_sd


boxplot(expr, las = 2, cex.axis = 0.7)


# ==========================================================
# Quantiles
# ==========================================================

quantile(expr[,1])


# ==========================================================
# Genes with low expression
# ==========================================================

lowest_expression <- counts %>%
  arrange(MeanExpression)

head(lowest_expression, 20)


# ==========================================================
# Genes with high expression
# ==========================================================

highest_expression <- counts %>%
  arrange(desc(MeanExpression))

head(highest_expression, 20)


# ==========================================================
# Saving results
# ==========================================================

write.csv(
  highest_expression,
  "results/highest_expression_genes.csv",
  row.names = FALSE
)

write.csv(
  lowest_expression,
  "results/lowest_expression_genes.csv",
  row.names = FALSE
)


# ==========================================================
# Lesson Summary
# ==========================================================

# Today you learned:
#
# ✓ How to explore an RNA-seq matrix
# ✓ How to summarize expression values
# ✓ How to visualize distributions
# ✓ How to identify highly expressed genes
# ✓ How to identify lowly expressed genes
# ✓ How to inspect sample variability
#
# ==========================================================


# ==========================================================
# Thinking Like a Bioinformatician
# ==========================================================

# Exploratory Data Analysis (EDA)
# is the first step of every RNA-seq analysis.
#
# Before asking:
#
# "Is there a significant difference?"
#
# Ask:
#
# "What do my data look like?"
#
# Good exploratory analysis prevents
# mistakes in downstream analyses.
#
# ==========================================================


# ==========================================================
# Challenge
# ==========================================================

# 1. Which sample has the highest average expression?
#
# 2. Which sample has the lowest average expression?
#
# 3. How many genes have MeanExpression < 1?
#
# 4. How many genes have MeanExpression > 100?
#
# 5. Which gene is the most highly expressed?
#
# 6. Which gene has the lowest expression?
#
# ==========================================================