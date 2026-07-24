# ==========================================================
#
# RNA-Seq Case Study
# Dataset: GSE245354
#
# Lesson 02
# Importing and Exploring the Expression Matrix
#
# Script:
# 02_import_expression_matrix.R
#
# Author: Jean Resende
#
# ==========================================================


# ==========================================================
# Lesson Objectives
# ==========================================================

# In this lesson you will learn:
#
# • How to locate the downloaded files
# • How to import an expression matrix into R
# • What is a data.frame
# • How to inspect an object
# • How to explore rows and columns
# • How to identify genes and samples
#
# ==========================================================


# ==========================================================
# Loading required packages
# ==========================================================

# No additional packages are required for this lesson.


# ==========================================================
# Locating the downloaded files
# ==========================================================

# Before importing the data, let's verify which files were
# downloaded in the previous lesson.

list.files(
  path = "data/raw/GSE245354",
  full.names = TRUE
)


# ==========================================================
# Defining the file path
# ==========================================================

# Instead of writing the complete file path several times,
# we store it inside an object.

expression_file <-
  "data/raw/GSE245354/GSE245354_Meegan_Normalized_counts.txt.gz"


# Display the file path

expression_file


# ==========================================================
# Importing the expression matrix
# ==========================================================

# The read.delim() function imports tab-delimited text files.
#
# The file is compressed (.gz), but R can read it directly
# without requiring manual decompression.

counts <- read.delim(expression_file)


# ==========================================================
# Inspecting the imported object
# ==========================================================

# What kind of object did R create?

class(counts)

# What is the internal structure?

str(counts)


# ==========================================================
# Matrix dimensions
# ==========================================================

# Number of rows and columns

dim(counts)

# Number of rows

nrow(counts)

# Number of columns

ncol(counts)


# ==========================================================
# Column names
# ==========================================================

# Display the names of every column.

names(counts)


# ==========================================================
# Previewing the dataset
# ==========================================================

# Display the first six rows.

head(counts)

# Display the last six rows.

tail(counts)


# ==========================================================
# Opening the data viewer
# ==========================================================

# Open the expression matrix using the RStudio Viewer.

View(counts)


# ==========================================================
# Exploring the dataset
# ==========================================================

# Which columns are available?

colnames(counts)

# Display the first column.

head(counts[,1])

# Display the second column.

head(counts[,2])

# Display the first six rows and first six columns.

counts[1:6, 1:6]


# ==========================================================
# Basic descriptive information
# ==========================================================

# Obtain a statistical summary of all numeric columns.

summary(counts)


# ==========================================================
# Lesson Summary
# ==========================================================

# In this lesson you learned:
#
# ✓ How to locate downloaded files
# ✓ How to import a tab-delimited file
# ✓ What is a data.frame
# ✓ How to inspect an object
# ✓ How to determine the dimensions of a dataset
# ✓ How to explore rows and columns
# ✓ How to preview an expression matrix


# ==========================================================
# Challenge
# ==========================================================

# 1. How many genes are present in the dataset?
#
# 2. How many samples are present?
#
# 3. What is the name of the first sample?
#
# 4. What is the name of the last sample?
#
# 5. Which column contains the gene identifiers?
#
# 6. Open the dataset using View() and explore it manually.
#
# 7. Can you identify which columns correspond to biological
#    samples?
#
# ==========================================================