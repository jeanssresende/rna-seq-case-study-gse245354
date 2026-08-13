# ==========================================================
#
# RNA-Seq Case Study
# Dataset: GSE245354
#
# Lesson 07
# Integrating Gene Expression with Sample Metadata
#
# Script:
# 07_integrating_expression_and_metadata.R
#
# Author: Jean Resende
#
# ==========================================================


# ==========================================================
# Lesson Objectives
# ==========================================================

# In this lesson you will learn:
#
# • Wide vs Long data format
# • How to reshape an expression matrix
# • How to integrate metadata
# • How to validate joins
# • How to prepare data for downstream analyses
#
# ==========================================================


# ==========================================================
# Loading packages
# ==========================================================

library(dplyr)
library(tidyr)
library(stringr)


# ==========================================================
# Loading expression matrix
# ==========================================================
counts <- read.delim("data/raw/GSE245354/GSE245354_Meegan_Normalized_counts.txt.gz")

# ==========================================================
# Loading metadata
# ==========================================================
metadata <- read.csv("data/metadata/sample_metadata.csv")

# ==========================================================
# Understanding the expression matrix
# ==========================================================
dim(counts)
head(counts)
tail(counts)

# ==========================================================
# Wide format
# ==========================================================

# Each sample occupies one column.
#
# This is called a WIDE table.
#
# Example:
#
# Gene   Sample1 Sample2 Sample3
# TP53      10      12      15
# ACTB     100      98     101
#
# ==========================================================


# ==========================================================
# Converting the entire matrix to LONG format
# ==========================================================

expression_long <- counts %>%
  pivot_longer(
    cols = 9:ncol(counts),
    names_to = "Sample",
    values_to = "Expression"
  )

expression_long
dim(expression_long)
head(expression_long)


# ==========================================================
# Why is this table so large?
# ==========================================================

# Number of genes
nrow(counts)

# Number of samples
ncol(counts) - 8

# Number of observations
nrow(expression_long)

# Expected:
# 25,506 genes × 24 samples = 612,144 observations


# ==========================================================
# Preparing metadata for integration
# ==========================================================

# The expression matrix uses sample names such as:
# X5958.NP.1_S1
#
# The GEO metadata contain:
# [5958-NP-1]
#
# We need to build a common identifier.

metadata <- metadata %>%
  mutate(
    MatrixSample = str_extract(Description, "5958-NP-\\d+"),
    MatrixSample = str_replace_all(MatrixSample, "-", "."),
    MatrixSample = paste0("X", MatrixSample, "_S1")
  )


metadata %>%
  select(Sample, MatrixSample, Treatment, Time, Donor) %>%
  head()

# ==========================================================
# Selecting one gene
# ==========================================================

gene_name = "MALAT1"

gene_expression <- expression_long %>%
  dplyr::filter(Gene.Symbol == gene_name)

gene_expression

dim(gene_expression)

# Expected: 24 rows


# ==========================================================
# Integrating metadata
# ==========================================================

gene_metadata <- gene_expression %>%
  left_join(
    metadata %>%
      select(MatrixSample, Treatment, Time, Donor),
    by = c("Sample" = "MatrixSample")
  )

gene_metadata

# ==========================================================
# Validating the join
# ==========================================================

sum(is.na(gene_metadata$Treatment))
sum(is.na(gene_metadata$Time))
sum(is.na(gene_metadata$Donor))

if (any(is.na(gene_metadata$Treatment))) {
  stop("Metadata integration failed: some samples were not matched.")
}


# ==========================================================
# Inspecting the integrated table
# ==========================================================

head(gene_metadata)
str(gene_metadata)
summary(gene_metadata)


# ==========================================================
# Expression by treatment
# ==========================================================

gene_metadata %>%
  group_by(Treatment) %>%
  summarise(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    N = n(),
    .groups = "drop"
  )


# ==========================================================
# Expression by collection time
# ==========================================================

gene_metadata %>%
  group_by(Time) %>%
  summarise(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    N = n(),
    .groups = "drop"
  )


# ==========================================================
# Expression by donor
# ==========================================================

gene_metadata %>%
  group_by(Donor) %>%
  summarise(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    N = n(),
    .groups = "drop"
  )


# ==========================================================
# Expression by Treatment and Time
# ==========================================================

gene_metadata %>%
  group_by(Treatment, Time) %>%
  summarise(
    Mean = mean(Expression),
    Median = median(Expression),
    SD = sd(Expression),
    N = n(),
    .groups = "drop"
  )


# ==========================================================
# Checking the experimental design
# ==========================================================

count(gene_metadata, Treatment)

count(gene_metadata, Time)

count(gene_metadata, Donor)


table(
  gene_metadata$Treatment,
  gene_metadata$Time
)


# ==========================================================
# Saving results
# ==========================================================

write.csv(
  gene_metadata,
  paste0("results/", tolower(gene_name), "_expression_metadata.csv"),
  row.names = FALSE)


# ==========================================================
# Lesson Summary
# ==========================================================

# Today you learned:
#
# ✓ Wide vs Long tables
# ✓ pivot_longer()
# ✓ left_join()
# ✓ Validation of joins
# ✓ group_by()
# ✓ summarise()
# ✓ Integration of metadata and expression
#
# ==========================================================


# ==========================================================
# Thinking Like a Bioinformatician
# ==========================================================

# RNA-seq analyses do not begin with
# statistical tests.
#
# They begin with data organization.
#
# A well-organized dataset makes every
# downstream analysis easier:
#
# • Boxplots
# • Violin plots
# • PCA
# • Heatmaps
# • Statistical tests
#
# The first step is always integrating
# expression data with metadata.
#
# Always validate your joins.
#
# A script that runs without errors is not
# necessarily a correct analysis.
#
# ==========================================================


# ==========================================================
# Challenge
# ==========================================================

# 1. Choose another gene.
#
# 2. Build the integrated table.
#
# 3. Compare expression among treatments.
#
# 4. Compare expression among collection times.
#
# 5. Compare expression among biological donors.
#
# 6. Verify that all joins returned zero NAs.
#
# ==========================================================