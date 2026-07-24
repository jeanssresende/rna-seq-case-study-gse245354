# ==========================================================
#
# RNA-Seq Case Study
# Dataset: GSE245354
#
# Lesson 01
# Downloading RNA-Seq Data from GEO
#
# Author: Jean Resende
#
# Repository:
# https://github.com/<your-user>/rna-seq-case-study-gse245354
#
# ==========================================================


# ==========================================================
# Lesson Objectives
# ==========================================================

# In this lesson you will learn:
#
# • What is the Gene Expression Omnibus (GEO)
# • How GEO organizes transcriptomic experiments
# • How to install Bioconductor packages
# • How to create objects in R
# • How to create a reproducible project structure
# • How to download public RNA-Seq datasets
# • How to inspect downloaded files


# ==========================================================
# About this dataset
# ==========================================================

# Throughout this course we will work with a real RNA-Seq
# dataset deposited in the NCBI Gene Expression Omnibus (GEO).
#
# GEO is one of the largest public repositories for
# functional genomics data and allows researchers to share
# sequencing experiments with the scientific community.
#
# Every experiment receives a unique accession number.
#
# Example:
#
# GSE245354
#
# where:
#
# GSE = GEO Series
#
# Each GEO Series contains one or more biological samples,
# identified by accession numbers beginning with GSM.


# ==========================================================
# Installing required packages
# ==========================================================

# BiocManager is the official package manager for
# Bioconductor.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


# GEOquery provides functions to download GEO datasets.

if (!requireNamespace("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")


# ==========================================================
# Loading libraries
# ==========================================================

library(GEOquery)


# ==========================================================
# Creating objects
# ==========================================================

# Objects store information in memory.
#
# During the analysis we will use the object below instead
# of repeatedly typing the accession number.

geo_accession <- "GSE245354"


# Display the object

geo_accession


# What type of object is this?

class(geo_accession)

typeof(geo_accession)

length(geo_accession)


# ==========================================================
# Creating the project structure
# ==========================================================

# Good bioinformatics projects separate raw data from
# processed data.
#
# The following commands create the folders used
# throughout this course.

dir.create("data", showWarnings = FALSE)

dir.create(
  "data/raw",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "data/metadata",
  showWarnings = FALSE
)

dir.create(
  "data/processed",
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)

dir.create(
  "results",
  showWarnings = FALSE
)


# ==========================================================
# Downloading supplementary files
# ==========================================================

# GEOquery automatically downloads all supplementary files
# associated with the selected GEO Series.

getGEOSuppFiles(
  GEO = geo_accession,
  baseDir = "data/raw"
)


# ==========================================================
# Inspecting downloaded files
# ==========================================================

# List every downloaded file.

list.files(
  path = "data/raw",
  recursive = TRUE
)


# ==========================================================
# Lesson Summary
# ==========================================================

# In this lesson you learned:
#
# ✓ What GEO is
# ✓ What a GEO Series is
# ✓ How to install Bioconductor packages
# ✓ How to create objects
# ✓ How to inspect objects
# ✓ How to organize a bioinformatics project
# ✓ How to download public RNA-Seq data


# ==========================================================
# Challenge
# ==========================================================

# 1. Visit the GEO website.
#
# https://www.ncbi.nlm.nih.gov/geo/
#
# 2. Search for another RNA-Seq dataset.
#
# 3. Replace the accession number used in this script.
#
# 4. Download the supplementary files.
#
# 5. Compare the downloaded files with GSE245354.
#
# Which experiment contains more samples?