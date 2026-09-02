# RNA-Seq Case Study – GSE245354

[![R](https://img.shields.io/badge/R-4.5+-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.22-019733?logo=bioconductor&logoColor=white)](https://bioconductor.org/)
[![RNA-Seq](https://img.shields.io/badge/RNA--Seq-Transcriptomics-blueviolet)]()
[![GEO](https://img.shields.io/badge/GEO-GSE245354-success)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245354)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-In%20Development-orange)]()

---

## Overview

This repository contains a practical **RNA-Seq case study** developed to teach **R programming, data analysis, and bioinformatics** using a real public transcriptomic dataset from the **NCBI Gene Expression Omnibus (GEO)**.

The project uses **GSE245354** as the main RNA-Seq case study and follows a progressive, lesson-based workflow. The objective is not only to teach R syntax, but to show how computational methods are applied to real biological data.

The analysis progresses from downloading and inspecting the dataset to expression exploration, metadata construction, data integration, PCA, heatmap analysis, validation of genes of interest, visualization, and statistical comparison.

The repository is also being expanded to include additional transcriptomic case studies, including **GSE25014**, a microarray dataset. The original GSE245354 analysis remains the main workflow and is not being replaced.

---

## Main Case Study: GSE245354

**Accession:** GSE245354  
**Data type:** RNA-Seq expression data  
**Source:** NCBI Gene Expression Omnibus (GEO)

The expression matrix used in the lessons contains **normalized Counts Per Million (CPM)** values rather than raw sequencing counts. This distinction is explicitly investigated during the course because the choice of statistical method depends on the type of expression data available.

The dataset contains **24 samples**, organized according to treatment, collection time, and biological donor.

### Experimental design

- **4 treatments:** Ctrl, Hb2+, Hb3+, Hemin
- **2 collection times:** 30 min, 1.5 hr
- **3 biological donors**
- **24 samples total**
- Balanced factorial design: **4 × 2 × 3 = 24**

This design provides an excellent example for teaching the relationship between an **expression matrix** and its corresponding **sample metadata**.

---

## Biological and Computational Questions

Throughout the case study, students progressively answer questions such as:

- What is GEO and how are transcriptomic experiments organized?
- What is an expression matrix?
- Are the values raw counts or normalized measurements?
- What does CPM mean?
- How can we construct metadata from GEO information?
- How do we connect expression data with experimental variables?
- How similar are the samples globally?
- Does treatment, time, or donor explain sample variability?
- Which genes are relevant to the biological question?
- How should multiple comparisons be handled?
- How can exploratory results be connected to statistical evidence?

---

## Course Workflow

```text
                    GSE245354
                        │
                        ▼
              Download data from GEO
                        │
                        ▼
              Import expression matrix
                        │
                        ▼
             Investigate expression data
                        │
                        ▼
              Build sample metadata
                        │
                        ▼
          Data manipulation with dplyr
                        │
                        ▼
          Exploratory data analysis (EDA)
                        │
                        ▼
          Integrate expression + metadata
                        │
                        ▼
                 PCA + Heatmap
                        │
                        ▼
             Validate target genes
                        │
                        ▼
            Visualize target genes
                        │
                        ▼
              Statistical analysis
                        │
                        ▼
             Biological interpretation
```

---

## Lessons

The repository currently contains the following lessons and scripts:

| Lesson | Script | Main topic |
|:------:|--------|------------|
| 01 | `01_download_data.R` | Downloading RNA-Seq data from GEO and organizing the project |
| 02 | `02_import_expression_matrix.R` | Importing and exploring the expression matrix |
| 03 | `03_investigating_rnaseq_data.R` | Investigating normalized expression values and understanding CPM |
| 04 | `04_data_manipulation_dplyr.R` | Data manipulation and transformation with dplyr |
| 05 | `05_building_metadata.R` | Building sample metadata from GEO information |
| 06 | `06_exploratory_data_analysis.R` | Exploratory data analysis and expression distributions |
| 07 | `07_integrating_expression_and_metadata.R` | Integrating expression data with metadata and validating joins |
| 08 | `08_pca_and_heatmap_analysis.R` | PCA, highly variable genes, heatmaps, and sample clustering |
| 09 | `09_genes_interesse_validacao.R` | Validation and extraction of 20 genes of interest |
| 10 | `10_genes_interesse_visualizacao.R` | Exploratory visualization of target genes |
| 11 | `11_genes_interesse_estatistica.R` | Statistical comparisons, paired tests, FDR, and effect exploration |

---

## From Expression Matrix to Biological Question

One of the central ideas of this project is that bioinformatics analysis is not simply about producing plots.

The workflow connects three essential components:

### 1. Expression matrix

```text
Genes × Samples
```

This table tells us **how much each gene is expressed** in each sample.

### 2. Sample metadata

```text
Samples × Experimental variables
```

This table tells us **what each sample represents**, including treatment, time, and donor.

### 3. Statistical analysis

The statistical analysis connects the biological question to the experimental design while accounting for important sources of variability.

This separation between expression data, metadata, and statistical modeling is a central principle of reproducible transcriptomic analysis.

---

## Genes of Interest

In Lessons 09–11, the project focuses on a predefined panel of **20 genes of biological interest**:

```text
CDH5
TJP1
CTTN
PTK2
HMOX1
NFE2L2
GPX4
SLC7A11
ACSL4
TFRC
FTH1
NCOA4
AIFM2
SOD2
NQO1
AKT1
MAPK1
MAPK3
BMPR2
TGFB1
```

The workflow teaches students how to:

- verify whether genes are present in the expression matrix;
- detect duplicated gene symbols;
- verify how many rows represent each gene;
- extract a target-gene expression matrix;
- integrate the genes with sample metadata;
- generate exploratory plots;
- compare expression between experimental conditions;
- apply multiple-testing correction;
- interpret statistical results in the context of the biological question.

---

## Statistical Analysis

The current statistical lesson uses the experimental structure of the dataset to introduce **paired comparisons** between treatments across the three biological donors.

The workflow includes:

- log2 transformation of CPM values;
- gene-wise Z-score calculation for visualization;
- paired t-tests;
- pairwise treatment comparisons;
- Benjamini–Hochberg FDR correction;
- effect-size exploration;
- volcano plots;
- expression plots with statistical annotations.

The donor is treated as an important source of biological variability because the same donors are represented across the experimental conditions.

The statistical workflow is presented as a teaching case and should be critically evaluated before being used as a final publication-grade model. Future lessons may extend this section to more formal linear-model or mixed-model approaches.

---

## Skills Covered

- R Programming
- RStudio
- Bioconductor
- GEOquery
- Data Wrangling
- dplyr
- tidyr
- stringr
- RNA-Seq
- Transcriptomics
- CPM and normalized expression
- Metadata construction
- Regular expressions
- Wide and long data formats
- Data integration and joins
- Exploratory Data Analysis
- Principal Component Analysis (PCA)
- Heatmaps
- Hierarchical clustering
- Gene-level exploration
- Statistical testing
- Multiple-testing correction
- Data visualization
- Reproducible Research
- Bioinformatics

---

## Repository Structure

```text
rna-seq-case-study-gse245354
│
├── README.md
├── LICENSE
├── rna-seq-case-study-gse245354.Rproj
│
├── data
│   ├── raw
│   ├── metadata
│   └── processed
│
├── scripts
│   ├── 01_download_data.R
│   ├── 02_import_expression_matrix.R
│   ├── 03_investigating_rnaseq_data.R
│   ├── 04_data_manipulation_dplyr.R
│   ├── 05_building_metadata.R
│   ├── 06_exploratory_data_analysis.R
│   ├── 07_integrating_expression_and_metadata.R
│   ├── 08_pca_and_heatmap_analysis.R
│   ├── 09_genes_interesse_validacao.R
│   ├── 10_genes_interesse_visualizacao.R
│   └── 11_genes_interesse_estatistica.R
│
├── figures
│
├── results
│
└── docs
```

Large raw and processed datasets are not intended to be version-controlled unless explicitly required for reproducibility.

---

## Getting Started

Clone the repository:

```bash
git clone https://github.com/jeanssresende/rna-seq-case-study-gse245354.git
```

Open the project in **RStudio** and run the scripts sequentially, beginning with:

```text
scripts/01_download_data.R
```

Each lesson builds on concepts and files generated in the previous lessons.

---

## Main R Packages

The project uses packages from CRAN and Bioconductor, including:

- GEOquery
- BiocManager
- tidyverse
- dplyr
- tidyr
- stringr
- ggplot2
- pheatmap
- RColorBrewer
- rstatix
- ggrepel

Additional packages may be introduced as the course evolves.

---

## Educational Philosophy

This repository was designed to bridge **R programming and bioinformatics** through a real-world transcriptomics case study.

The lessons follow a progressive approach:

```text
Learn the R concept
        ↓
Apply it to real biological data
        ↓
Inspect the result
        ↓
Ask a biological question
        ↓
Validate the analysis
        ↓
Interpret the result
```

The goal is to teach students to think like bioinformaticians rather than simply execute commands.

A central principle of the course is:

> **A script that runs without errors is not necessarily a correct analysis.**

Students are therefore encouraged to inspect objects, validate joins, investigate metadata, question assumptions, and understand the biological meaning of each analytical step.

---

## Future Case Studies

The repository will gradually expand beyond GSE245354 while preserving the analyses already developed for this case study.

### GSE25014 — Microarray

The next case study will introduce **GSE25014**, an expression-profiling microarray dataset. Unlike GSE245354, this dataset uses an Affymetrix platform and therefore introduces concepts that are specific to microarray analysis, including:

- probe-level expression;
- platform annotation;
- probe-to-gene mapping;
- multiple probes representing the same gene;
- probe summarization;
- construction of a gene-level expression matrix.

This second dataset will complement the RNA-Seq case study and allow students to compare two major transcriptomics technologies within the same educational repository.

**Important:** GSE25014 is an additional case study. It does **not** replace GSE245354.

---

## Learning Objectives

By completing the current GSE245354 case study, students should be able to:

- organize a reproducible bioinformatics project;
- obtain public transcriptomic data from GEO;
- understand the structure of an RNA-Seq expression matrix;
- distinguish raw counts from normalized expression values;
- understand CPM normalization;
- inspect GEO metadata;
- construct metadata tables programmatically;
- use regular expressions to extract experimental information;
- manipulate data with dplyr and tidyr;
- convert between wide and long data formats;
- integrate expression data with sample metadata;
- validate joins and detect missing matches;
- perform exploratory data analysis;
- evaluate sample relationships with PCA;
- construct and interpret heatmaps;
- identify and validate genes of interest;
- visualize expression patterns;
- perform paired statistical comparisons;
- apply multiple-testing correction;
- connect statistical results to biological questions.

---

## Author

**Jean Resende**  
PhD Student in Biotechnology  
Bioinformatics • Transcriptomics • Immunoinformatics

---

## License

This project is distributed under the MIT License.

---

## Citation

If you use this repository for teaching or research purposes, please cite the original **GSE245354** GEO record and its associated publication.

For additional case studies added to this repository, cite the corresponding GEO accession and original publication.
