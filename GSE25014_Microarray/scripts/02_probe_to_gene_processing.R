# ============================================================
# Transcriptomics Case Studies
# Case Study: GSE25014
# Technology: Microarray
# ============================================================
#
# Lesson 02: Probe-to-Gene Processing
#
# Dataset:
# GSE25014 - Heme Effects on Pulmonary Endothelial Cells
#
# Platform:
# Affymetrix Human Genome U133 Plus 2.0
#
# Objective:
#   Understand the relationship between microarray probes
#   and genes and investigate why multiple probes can
#   represent the same gene.
#
# Workflow:
#
#   Probe-level expression
#            ↓
#   Probe annotation
#            ↓
#   Probe-to-gene mapping
#            ↓
#   Multiple probes / gene
#            ↓
#   Gene-level expression
#
# IMPORTANT:
#   In this lesson we investigate and prepare the data.
#   The final probe summarization strategy will be
#   implemented in a subsequent lesson.
#
# ============================================================


# ============================================================
# Lesson Objectives
# ============================================================
#
# In this lesson you will learn:
#
# • How to load the data generated in Lesson 01
# • How to verify probe IDs between datasets
# • How to investigate probe annotation
# • How to identify probes without gene annotation
# • How to identify genes represented by multiple probes
# • How to investigate a specific gene
# • How to connect probe annotation with expression values
# • Why probe-level data need to be processed before
#   gene-level analysis
#
# ============================================================


# ============================================================
# 0. Load Required Packages
# ============================================================

library(tidyverse)


# ============================================================
# 1. Define Data Directory
# ============================================================

data_dir <- "GSE25014_Microarray/data/raw/GSE25014"


# ============================================================
# 2. Load Data from Lesson 01
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Loading GSE25014 data\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Expression matrix
# ------------------------------------------------------------
#
# Rows    → probes
# Columns → samples
#
# Expression values are RMA-normalized and log2 transformed.

expr <- read.csv(
  file.path(
    data_dir,
    "expr_matrix_rma.csv"
  ),
  row.names = 1,
  check.names = FALSE
)


# ------------------------------------------------------------
# Sample metadata
# ------------------------------------------------------------

metadata <- read.csv(
  file.path(
    data_dir,
    "metadata.csv"
  ),
  check.names = FALSE
)


# ------------------------------------------------------------
# Probe annotation
# ------------------------------------------------------------

annotations <- read.csv(
  file.path(
    data_dir,
    "probe_annotations.csv"
  ),
  check.names = FALSE
)


# ============================================================
# 3. Inspect Loaded Objects
# ============================================================

cat("\n")
cat("Expression matrix:\n")

print(
  class(expr)
)

print(
  dim(expr)
)


cat("\nMetadata:\n")

print(
  class(metadata)
)

print(
  dim(metadata)
)


cat("\nProbe annotation:\n")

print(
  class(annotations)
)

print(
  dim(annotations)
)


# ============================================================
# 4. Understand the Expression Matrix
# ============================================================
#
# A microarray expression matrix has the structure:
#
#
#              Samples
#           S1   S2   S3   S4
#
# Probe 1   ...  ...  ...  ...
# Probe 2   ...  ...  ...  ...
# Probe 3   ...  ...  ...  ...
#
#
# Each row represents a probe.
#
# The first question is:
#
# "Are the probe IDs in the expression matrix
#  also present in the annotation table?"
#
# ============================================================


probe_expr <- rownames(expr)

probe_annotation <- annotations$ProbeID


cat("\n")
cat("============================================================\n")
cat("Probe ID Check\n")
cat("============================================================\n")


cat(
  "Probes in expression matrix:",
  length(probe_expr),
  "\n"
)


cat(
  "Probes in annotation table:",
  length(probe_annotation),
  "\n"
)


# ============================================================
# 5. Compare Expression Probes with Annotation Probes
# ============================================================

matched_probes <- sum(
  probe_expr %in% probe_annotation
)


unmatched_probes <- sum(
  !probe_expr %in% probe_annotation
)


cat(
  "\nExpression probes with annotation:",
  matched_probes,
  "\n"
)


cat(
  "Expression probes without annotation:",
  unmatched_probes,
  "\n"
)


# Percentage of matched probes

percentage_matched <- (
  matched_probes /
    length(probe_expr)
) * 100


cat(
  "Percentage of probes with annotation:",
  round(
    percentage_matched,
    2
  ),
  "%\n"
)


# ============================================================
# 6. Identify Unmatched Probe IDs
# ============================================================
#
# Which probes are present in the expression matrix
# but are absent from the annotation table?
#
# ============================================================

unmatched_probe_ids <- probe_expr[
  !probe_expr %in% probe_annotation
]


cat("\n")
cat(
  "Number of unmatched probe IDs:",
  length(unmatched_probe_ids),
  "\n"
)


cat("\nFirst unmatched probes:\n")

print(
  head(
    unmatched_probe_ids,
    20
  )
)


# ============================================================
# 7. Inspect Probe Annotation
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Probe Annotation\n")
cat("============================================================\n")


cat("\nFirst probes:\n")

print(
  annotations %>%
    select(
      ProbeID,
      `Gene symbol`,
      `Gene title`
    ) %>%
    head(10)
)


# ------------------------------------------------------------
# Random sample of probes
# ------------------------------------------------------------
#
# A random sample helps us inspect different examples
# instead of looking only at the first rows.
#
# ============================================================

set.seed(123)


random_probes <- annotations %>%
  slice_sample(
    n = 10
  ) %>%
  select(
    ProbeID,
    `Gene symbol`,
    `Gene title`
  )


cat("\nRandom sample of probes:\n")

print(
  random_probes
)


# ============================================================
# 8. Probe Annotation Summary
# ============================================================
#
# Not every probe necessarily has a gene symbol.
#
# We separate:
#
#   1. Probes with gene symbol
#   2. Probes without gene symbol
#
# ============================================================


probes_with_gene <- annotations %>%
  filter(
    !is.na(`Gene symbol`),
    `Gene symbol` != ""
  )


probes_without_gene <- annotations %>%
  filter(
    is.na(`Gene symbol`) |
      `Gene symbol` == ""
  )


cat("\n")
cat("============================================================\n")
cat("Probe Annotation Summary\n")
cat("============================================================\n")


cat(
  "\nProbes with gene symbol:",
  nrow(probes_with_gene),
  "\n"
)


cat(
  "Probes without gene symbol:",
  nrow(probes_without_gene),
  "\n"
)


# Create summary table

annotation_summary <- tibble(
  
  Category = c(
    "With gene symbol",
    "Without gene symbol"
  ),
  
  Probes = c(
    nrow(probes_with_gene),
    nrow(probes_without_gene)
  )
  
)


cat("\nAnnotation summary:\n")

print(
  annotation_summary
)


# ============================================================
# 9. Calculate Number of Probes per Gene
# ============================================================
#
# This is one of the most important steps of the lesson.
#
# In microarray data:
#
#        Probe 1 ──┐
#        Probe 2 ──┼──→ Gene A
#        Probe 3 ──┘
#
#
# Therefore:
#
#        1 gene ≠ necessarily 1 probe
#
# We need to identify genes represented by multiple probes.
#
# ============================================================


probe_count_per_gene <- annotations %>%
  
  filter(
    !is.na(`Gene symbol`),
    `Gene symbol` != ""
  ) %>%
  
  count(
    `Gene symbol`,
    name = "Number_of_probes"
  ) %>%
  
  arrange(
    desc(Number_of_probes)
  )


cat("\n")
cat("============================================================\n")
cat("Number of Probes per Gene\n")
cat("============================================================\n")


print(
  head(
    probe_count_per_gene,
    30
  )
)


# ============================================================
# 10. Summary of Probe Multiplicity
# ============================================================
#
# How many genes have:
#
#   1 probe?
#   2 probes?
#   3 probes?
#   ...
#
# ============================================================


probe_multiplicity <- probe_count_per_gene %>%
  
  count(
    Number_of_probes,
    name = "Number_of_genes"
  ) %>%
  
  arrange(
    Number_of_probes
  )


cat("\n")
cat("Distribution of probes per gene:\n")

print(
  probe_multiplicity
)


# ============================================================
# 11. Select a Gene for Investigation
# ============================================================
#
# We will first investigate CDH5.
#
# CDH5 is represented by multiple probes in the
# Affymetrix platform.
#
# The objective is NOT yet to choose the best probe.
#
# The objective is to understand:
#
#   Gene
#      ↓
#   Multiple probes
#      ↓
#   Different expression measurements
#
# ============================================================


gene_of_interest <- "CDH5"


cat("\n")
cat("============================================================\n")
cat("Gene Investigation:", gene_of_interest, "\n")
cat("============================================================\n")


# Find probes associated with CDH5

cdh5_probes <- annotations %>%
  
  filter(
    `Gene symbol` == gene_of_interest
  ) %>%
  
  select(
    ProbeID,
    `Gene symbol`,
    `Gene title`
  )


cat("\nProbes representing", gene_of_interest, ":\n")

print(
  cdh5_probes
)


cat(
  "\nNumber of probes:",
  nrow(cdh5_probes),
  "\n"
)


# ============================================================
# 12. Extract CDH5 Expression
# ============================================================

cdh5_probe_ids <- cdh5_probes$ProbeID


cdh5_expression <- expr[
  
  rownames(expr) %in% cdh5_probe_ids,
  
  ,
  
  drop = FALSE
  
]


cat("\n")
cat("CDH5 expression matrix:\n")

print(
  cdh5_expression
)


# ============================================================
# 13. Combine Probe Annotation and Expression
# ============================================================
#
# Now we combine:
#
#   Probe ID
#   Gene symbol
#   Gene title
#   Expression values
#
# This creates a table that allows us to inspect
# the relationship between probes and their expression.
#
# ============================================================


cdh5_table <- cdh5_expression %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    var = "ProbeID"
  ) %>%
  
  left_join(
    cdh5_probes,
    by = "ProbeID"
  ) %>%
  
  select(
    ProbeID,
    `Gene symbol`,
    `Gene title`,
    everything()
  )


cat("\n")
cat("CDH5 probe-expression table:\n")

print(
  cdh5_table
)


# ============================================================
# 14. Compare CDH5 Probe Expression
# ============================================================
#
# For each probe, calculate:
#
#   Mean expression across samples
#   Median expression across samples
#   Standard deviation
#
# IMPORTANT:
#
# These statistics are only descriptive.
#
# We are NOT choosing a probe yet.
#
# ============================================================


cdh5_probe_summary <- cdh5_expression %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    var = "ProbeID"
  ) %>%
  
  mutate(
    
    Mean = rowMeans(
      select(
        .,
        -ProbeID
      ),
      na.rm = TRUE
    ),
    
    Median = apply(
      select(
        .,
        -ProbeID
      ),
      1,
      median,
      na.rm = TRUE
    ),
    
    SD = apply(
      select(
        .,
        -ProbeID
      ),
      1,
      sd,
      na.rm = TRUE
    )
    
  ) %>%
  
  select(
    ProbeID,
    Mean,
    Median,
    SD
  ) %>%
  
  arrange(
    desc(Mean)
  )


cat("\n")
cat("CDH5 probe summary:\n")

print(
  cdh5_probe_summary
)


# ============================================================
# 15. Investigate a Second Gene: CD44
# ============================================================
#
# CD44 provides another example of a gene represented
# by multiple probes.
#
# ============================================================


gene_of_interest <- "CD44"


cat("\n")
cat("============================================================\n")
cat("Gene Investigation:", gene_of_interest, "\n")
cat("============================================================\n")


cd44_probes <- annotations %>%
  
  filter(
    `Gene symbol` == gene_of_interest
  ) %>%
  
  select(
    ProbeID,
    `Gene symbol`,
    `Gene title`
  )


cat("\nProbes representing", gene_of_interest, ":\n")

print(
  cd44_probes
)


cat(
  "\nNumber of probes:",
  nrow(cd44_probes),
  "\n"
)


# ============================================================
# 16. Extract CD44 Expression
# ============================================================


cd44_probe_ids <- cd44_probes$ProbeID


cd44_expression <- expr[
  
  rownames(expr) %in% cd44_probe_ids,
  
  ,
  
  drop = FALSE
  
]


cat("\n")
cat("CD44 expression matrix:\n")

print(
  cd44_expression
)


# ============================================================
# 17. Summarize CD44 Probes
# ============================================================


cd44_probe_summary <- cd44_expression %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    var = "ProbeID"
  ) %>%
  
  mutate(
    
    Mean = rowMeans(
      select(
        .,
        -ProbeID
      ),
      na.rm = TRUE
    ),
    
    Median = apply(
      select(
        .,
        -ProbeID
      ),
      1,
      median,
      na.rm = TRUE
    ),
    
    SD = apply(
      select(
        .,
        -ProbeID
      ),
      1,
      sd,
      na.rm = TRUE
    )
    
  ) %>%
  
  select(
    ProbeID,
    Mean,
    Median,
    SD
  ) %>%
  
  arrange(
    desc(Mean)
  )


cat("\n")
cat("CD44 probe summary:\n")

print(
  cd44_probe_summary
)


# ============================================================
# 18. Key Concept
# ============================================================
#
# At this point we have observed that:
#
#
#             MICROARRAY
#                 │
#                 ↓
#              PROBES
#                 │
#          annotation
#                 ↓
#               GENES
#
#
# However, a gene may have multiple probes:
#
#
#       Probe A ──┐
#       Probe B ──┼──→ Gene X
#       Probe C ──┘
#
#
# Each probe can have a different expression profile.
#
# Therefore, we cannot simply assume:
#
#       one probe = one gene
#
#
# Before performing gene-level analyses, we need to
# define a strategy for summarizing multiple probes.
#
# ============================================================


# ============================================================
# 19. What We Will Do in the Next Lesson
# ============================================================
#
# In the next lesson we will investigate:
#
# • Which probes should be retained?
# • How should probes without gene annotation be handled?
# • How should duplicated gene symbols be treated?
# • What strategies can be used to summarize multiple probes?
#
# Possible strategies include:
#
#   - Mean
#   - Median
#   - Maximum expression
#   - Maximum variance
#   - Probe-selection methods
#
# The choice of strategy should be biologically and
# statistically justified.
#
# ============================================================


# ============================================================
# Lesson Summary
# ============================================================
#
# ✓ Loaded expression, metadata and annotation data
#
# ✓ Verified probe IDs
#
# ✓ Identified matched and unmatched probes
#
# ✓ Investigated probe annotation
#
# ✓ Quantified probes with and without gene symbols
#
# ✓ Counted probes per gene
#
# ✓ Identified genes represented by multiple probes
#
# ✓ Investigated CDH5
#
# ✓ Investigated CD44
#
# ✓ Connected probe annotation with expression values
#
# ✓ Calculated descriptive statistics for probes
#
# ✓ Identified the need for probe summarization
#
# ============================================================


# ============================================================
# Challenge
# ============================================================
#
# Answer the following questions:
#
# 1. How many probes are present in the expression matrix?
#
# 2. How many probes have a gene symbol?
#
# 3. How many probes do not have a gene symbol?
#
# 4. How many genes are represented by more than one probe?
#
# 5. What is the gene with the largest number of probes?
#
# 6. How many probes represent CDH5?
#
# 7. How many probes represent CD44?
#
# 8. Do the probes representing the same gene have
#    identical expression values?
#
# 9. Why can different probes representing the same gene
#    produce different expression profiles?
#
# 10. Why can we not simply keep the first probe for
#     each gene?
#
# ============================================================


# ============================================================
# End of Lesson 02
# ============================================================