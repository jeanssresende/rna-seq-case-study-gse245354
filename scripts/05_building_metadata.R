# ==========================================================
#
# RNA-Seq Case Study
# Dataset: GSE245354
#
# Lesson 05
# Building Sample Metadata from GEO
#
# Script:
# 05_building_metadata.R
#
# Author: Jean Resende
#
# ==========================================================


# ==========================================================
# Lesson Objectives
# ==========================================================

# In this lesson you will learn:
#
# • What metadata are
# • Why metadata are essential
# • How to download metadata from GEO
# • How to manipulate character strings
# • How to use regular expressions
# • How to build a metadata table automatically
#
# ==========================================================


# ==========================================================
# Installing required packages
# ==========================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!requireNamespace("stringr", quietly = TRUE))
  install.packages("stringr")

library(GEOquery)
library(dplyr)
library(stringr)


# ==========================================================
# Downloading GEO metadata
# ==========================================================

gse <- getGEO("GSE245354", GSEMatrix = TRUE)

metadata <- pData(gse[[1]])


# ==========================================================
# Exploring the metadata
# ==========================================================

class(metadata)
dim(metadata)
names(metadata)
head(metadata)
View(metadata)


# ==========================================================
# Sample titles
# ==========================================================

# The title contains useful experimental information.
metadata$title

# Display unique sample descriptions.
unique(metadata$title)


# ==========================================================
# Creating a metadata table
# ==========================================================

metadata_clean <- metadata %>%
  select(Sample = geo_accession,
         Description = title)

metadata_clean


# ==========================================================
# Extracting the treatment
# ==========================================================

# Sample description:
#
# Ctrl 30 min 1
#
# Hb2+ 30 min 1
#
# Hb3+ 30 min 1
#
# Hemin 30 min 1

metadata_clean <- metadata_clean %>%
  mutate(Treatment = str_extract(Description, "Ctrl|Hb2\\+|Hb3\\+|Hemin"))

metadata_clean


# ==========================================================
# Extracting the collection time
# ==========================================================

metadata_clean <- metadata_clean %>%
  mutate( Time = str_extract(Description, "30 min|1\\.5 hr"))

metadata_clean


# ==========================================================
# Extracting the biological replicate
# ==========================================================

# We want to capture the number located
# immediately after "min" or "hr".

metadata_clean <- metadata_clean %>%
  mutate(Donor = str_extract(Description, "(?<=min )\\d|(?<=hr )\\d"))

metadata_clean


# ==========================================================
# Converting variables into factors
# ==========================================================

metadata_clean <- metadata_clean %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Ctrl", "Hb2+", "Hb3+", "Hemin")),
         Time = factor(Time,
                       levels = c("30 min", "1.5 hr")),
         Donor = factor(Donor))


str(metadata_clean)


# ==========================================================
# Checking the experimental design
# ==========================================================

count(metadata_clean, Treatment)

count(metadata_clean, Time)

count(metadata_clean, Donor)


table(metadata_clean$Treatment, metadata_clean$Time)

table(metadata_clean$Treatment, metadata_clean$Donor)


# ==========================================================
# Saving the metadata
# ==========================================================

write.csv(metadata_clean,
          file = "data/metadata/sample_metadata.csv",
          row.names = FALSE)


# ==========================================================
# Lesson Summary
# ==========================================================

# Today you learned:
#
# ✓ How to download metadata from GEO
#
# ✓ How to inspect sample information
#
# ✓ How to extract information from text
#
# ✓ How to use regular expressions
#
# ✓ How to build a metadata table
#
# ✓ How to identify the experimental design
#
# ✓ How to save metadata for future analyses


# ==========================================================
# Thinking Like a Bioinformatician
# ==========================================================

# A bioinformatics project is built upon two tables:
#
# 1) Expression Matrix
#
# Genes × Samples
#
#
# 2) Metadata
#
# Samples × Experimental Variables
#
#
# The expression matrix tells us
# HOW MUCH each gene is expressed.
#
# The metadata tell us
# WHAT each sample represents.
#
# Without metadata,
# biological interpretation is impossible.
#
#
# Never type manually information that can
# be extracted automatically.
#
# Automation makes analyses reproducible.
#
# ==========================================================


# ==========================================================
# Challenge
# ==========================================================

# 1. How many samples are there?
#
# 2. How many treatments?
#
# 3. How many biological replicates?
#
# 4. Which treatment has the most samples?
#
# 5. Is the experimental design balanced?
#
# 6. Open sample_metadata.csv.
#
# 7. Compare this table with the original GEO metadata.
#
# ==========================================================