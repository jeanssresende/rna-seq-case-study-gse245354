# ==========================================================
#
# RNA-Seq Case Study
# Dataset: GSE245354
#
# Lesson 04
# Manipulating Gene Expression Data with dplyr
#
# Script:
# 04_data_manipulation_dplyr.R
#
# Author: Jean Resende
#
# ==========================================================


# ==========================================================
# Lesson Objectives
# ==========================================================

# In this lesson you will learn:
#
# • What is the tidyverse
# • What is dplyr
# • How to use the pipe operator (%>%)
# • How to select columns
# • How to filter rows
# • How to rename variables
# • How to reorder columns
# • How to sort a dataset
#
# ==========================================================


# ==========================================================
# Installing and loading dplyr
# ==========================================================

# Install dplyr only if necessary.

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(dplyr)


# ==========================================================
# Importing the expression matrix
# ==========================================================

expression_file <-
  "data/raw/GSE245354/GSE245354_Meegan_Normalized_counts.txt.gz"

counts <- read.delim(expression_file)


# ==========================================================
# Understanding the Pipe Operator
# ==========================================================

# The pipe (%>%) passes the result of one operation
# directly to the next.
#
# Instead of writing:
#
# function2(function1(data))
#
# We can write:
#
# data %>%
#   function1() %>%
#   function2()
#
# This makes the code easier to read.


# ==========================================================
# Selecting Columns
# ==========================================================

# Which columns are available?
names(counts)


# Select only the gene symbol.
counts %>%
  select(Gene.Symbol)


# Select annotation columns.
counts %>%
  select(
    Chromosome,
    Gene.Symbol,
    gene_biotype
  )


# Select annotation plus the first three samples.
counts %>%
  select(
    Chromosome,
    Gene.Symbol,
    gene_biotype,
    X5958.NP.18_S1,
    X5958.NP.15_S1,
    X5958.NP.17_S1
  )


# ==========================================================
# Renaming Columns
# ==========================================================
counts %>%
  rename(
    Gene = Gene.Symbol,
    Biotype = gene_biotype
  )


# ==========================================================
# Reordering Columns
# ==========================================================
counts %>%
  relocate(
    gene_biotype,
    .after = Gene.Symbol
  )


# ==========================================================
# Filtering Rows
# ==========================================================

# Keep only protein coding genes.
counts %>%
  filter(
    gene_biotype == "protein_coding"
  )


# Genes located on chromosome 1.
counts %>%
  filter(
    Chromosome == "1"
  )


# Protein coding genes from chromosome 1.
counts %>%
  filter(
    Chromosome == "1",
    gene_biotype == "protein_coding"
  )


# ==========================================================
# Sorting the Dataset
# ==========================================================

# Lowest expression.
counts %>%
  arrange(
    X5958.NP.18_S1
  )


# Highest expression.
counts %>%
  arrange(
    desc(X5958.NP.18_S1)
  )


# ==========================================================
# Combining Multiple Verbs
# ==========================================================

counts %>%
  filter(
    gene_biotype == "protein_coding"
  ) %>%
  select(
    Gene.Symbol,
    X5958.NP.18_S1
  ) %>%
  arrange(
    desc(X5958.NP.18_S1)
  )


# ==========================================================
# Saving the Result
# ==========================================================

top_genes <-
  counts %>%
  filter(
    gene_biotype == "protein_coding"
  ) %>%
  select(
    Gene.Symbol,
    X5958.NP.18_S1
  ) %>%
  arrange(
    desc(X5958.NP.18_S1)
  )


head(top_genes)


# ==========================================================
# Lesson Summary
# ==========================================================

# Today you learned:
#
# ✓ What is dplyr
# ✓ How the pipe operator works
# ✓ How to select variables
# ✓ How to rename columns
# ✓ How to reorder columns
# ✓ How to filter observations
# ✓ How to sort a dataset
# ✓ How to combine multiple dplyr verbs


# ==========================================================
# Thinking Like a Bioinformatician
# ==========================================================

# Bioinformatics is not about memorizing functions.
#
# It is about asking biological questions.
#
# Every dplyr verb answers one type of question:
#
# select()   -> Which variables do I need?
#
# filter()   -> Which genes satisfy my criteria?
#
# arrange()  -> Which genes have the highest expression?
#
# rename()   -> Can I make the dataset easier to read?
#
# relocate() -> Can I improve the organization of my table?
#
# Learning dplyr means learning how to interrogate
# biological data.


# ==========================================================
# Challenge
# ==========================================================

# 1. Select only Gene.Symbol and gene_biotype.
#
# 2. Display only genes located on chromosome X.
#
# 3. Display only lncRNA genes.
#
# 4. Rename Gene.Symbol to Gene.
#
# 5. Sort chromosome 1 genes by decreasing expression
#    in sample GSM7841265.
#
# 6. Create an object called chromosome1_genes.
#
# 7. Display the first 20 genes.
#
# ==========================================================