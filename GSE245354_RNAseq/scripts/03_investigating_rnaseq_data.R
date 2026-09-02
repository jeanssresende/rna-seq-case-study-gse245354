# ==========================================================
#
# RNA-Seq Case Study
# Dataset: GSE245354
#
# Lesson 03
# Investigating RNA-Seq Expression Data
#
# Script:
# 03_investigating_rnaseq_data.R
#
# Author: Jean Resende
#
# ==========================================================


# ==========================================================
# Lesson Objectives
# ==========================================================

# In this lesson you will learn:
#
# • How to investigate an RNA-Seq expression matrix
# • How to recognize normalized expression values
# • What Counts Per Million (CPM) means
# • How to use GEO metadata
# • Why preprocessing information is important
# • Why we should understand the data before performing
#   statistical analyses
#
# ==========================================================


# ==========================================================
# Importing the expression matrix
# ==========================================================

expression_file <-
  "data/raw/GSE245354/GSE245354_Meegan_Normalized_counts.txt.gz"

counts <- read.delim(expression_file)


# ==========================================================
# Investigation 1
# What kind of numbers are these?
# ==========================================================

# Let's inspect the first rows.

head(counts)

# ----------------------------------------------------------
# Question
#
# Do these values represent sequencing reads?
#
# Example:
#
# 5.734850
#
# Could we have 5.734850 sequencing reads?
#
# Think before continuing.
#
# ----------------------------------------------------------


# ==========================================================
# Investigation 2
# Are these raw counts?
# ==========================================================

# Raw counts should be integers.
#
# Let's verify whether all expression values are integers.

expression_matrix <- counts[, 9:ncol(counts)]

all(expression_matrix %% 1 == 0)

# ----------------------------------------------------------
# Expected answer:
#
# FALSE
#
# Therefore, these are not raw counts.
# ----------------------------------------------------------


# ==========================================================
# Investigation 3
# Minimum and maximum values
# ==========================================================

range(expression_matrix)

summary(expression_matrix)

# ----------------------------------------------------------
# Questions
#
# Why is the minimum value 0.0001?
#
# Why are there decimal values?
#
# ----------------------------------------------------------


# ==========================================================
# Investigation 4
# Are there negative values?
# ==========================================================

any(expression_matrix < 0)

# Negative values usually indicate transformed data.
#
# In this dataset, all values are positive.


# ==========================================================
# Investigation 5
# Average expression per sample
# ==========================================================

colMeans(expression_matrix)

# Observe that all samples have very similar means.
#
# Why might this happen?


# ==========================================================
# Investigation 6
# Visualizing the distribution
# ==========================================================

hist(
  unlist(expression_matrix),
  breaks = 100,
  main = "Distribution of Expression Values",
  xlab = "Expression"
)


# ==========================================================
# Investigation 7
# Boxplot
# ==========================================================

boxplot(
  expression_matrix,
  las = 2,
  main = "Expression Distribution Across Samples"
)


# ==========================================================
# Looking outside R
# ==========================================================

# Bioinformatics does not end inside R.
#
# We must investigate the metadata available in GEO.
#
# Open the GEO page and inspect the "Data processing"
# section of one sample.


# ==========================================================
# GEO Metadata
# ==========================================================

# Data processing:
#
# TopHat2 was used to align reads.
#
# Reads were quantified using the Partek E/M annotation model.
#
# Reads were normalized to counts per million (CPM)
# with 0.0001 added to the result.
#
# Supplementary files include CPM values for each sample.
#
# ----------------------------------------------------------


# ==========================================================
# Conclusion
# ==========================================================

# The expression matrix contains:
#
# ✓ Counts Per Million (CPM)
#
# It does NOT contain raw read counts.
#
# A pseudocount of 0.0001 was added to every value.


# ==========================================================
# What is CPM?
# ==========================================================

# CPM = Counts Per Million
#
#                   Gene Counts
# CPM = -------------------------------- × 1,000,000
#        Total Reads in the Sample
#
# CPM corrects differences in sequencing depth between
# libraries, allowing expression values from different
# samples to be compared more fairly.


# ==========================================================
# Why was 0.0001 added?
# ==========================================================

# Possible reasons:
#
# • Avoid zero values
# • Prevent mathematical problems
# • Facilitate logarithmic transformations
# • Improve downstream analyses


# ==========================================================
# Can we use DESeq2?
# ==========================================================

# Important:
#
# DESeq2 requires RAW COUNTS.
#
# This matrix contains normalized CPM values.
#
# Therefore:
#
# This matrix should NOT be used directly in DESeq2.
#
# If differential expression analysis is required,
# raw counts must be obtained.


# ==========================================================
# Lesson Summary
# ==========================================================

# Today you learned:
#
# ✓ How to investigate an RNA-Seq dataset
# ✓ How to recognize normalized data
# ✓ What CPM means
# ✓ Why metadata are essential
# ✓ Why we should never trust only the file name
# ✓ How GEO documents preprocessing steps
# ✓ Why raw counts and normalized counts are different


# ==========================================================
# Challenge
# ==========================================================

# 1. What evidence indicates that this matrix is not composed
#    of raw counts?
#
# 2. Which software was used for read alignment?
#
# 3. Which normalization method was applied?
#
# 4. Why was 0.0001 added to the expression values?
#
# 5. What does CPM stand for?
#
# 6. Can this matrix be used directly with DESeq2?
#
# 7. Explain your answer.
#
# ==========================================================