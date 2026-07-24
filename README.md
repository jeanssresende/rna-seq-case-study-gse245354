# RNA-Seq Case Study – GSE245354

[![R](https://img.shields.io/badge/R-4.5+-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.22-019733?logo=bioconductor&logoColor=white)](https://bioconductor.org/)
[![Linux](https://img.shields.io/badge/Linux-Ubuntu-E95420?logo=ubuntu&logoColor=white)](https://ubuntu.com/)
[![RNA-Seq](https://img.shields.io/badge/RNA--Seq-Transcriptomics-blueviolet)]()
[![GEO](https://img.shields.io/badge/GEO-GSE245354-success)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245354)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-In%20Development-orange)]()

---

## Overview

This repository contains a practical **RNA-Seq case study** developed to teach **R programming for transcriptomic data analysis** using a real public dataset from the **NCBI Gene Expression Omnibus (GEO)**.

Rather than presenting isolated R commands, this project follows the workflow adopted in real bioinformatics research, allowing students to learn programming while performing a complete RNA-Seq exploratory analysis.

The lessons are organized as a sequence of scripts that progressively introduce data import, manipulation, exploratory analysis, statistical inference, visualization, and interpretation of transcriptomic data.

---

## Getting Started

Clone this repository:

```bash
git clone https://github.com/jeanssresende/rna-seq-case-study-gse245354.git
```

Open the project in RStudio and execute:

```r
scripts/01_download_data.R
```

The script will automatically:

- install required packages (if necessary);
- create the project directory structure;
- download the GSE245354 dataset into `data/raw/`.

---

## Skills Covered

- R Programming
- Data Wrangling
- Tidyverse
- RNA-Seq
- Transcriptomics
- Statistics
- Data Visualization
- Reproducible Research
- Bioinformatics
- Computational Biology

---

## Dataset

**Accession:** GSE245354

**Source:** NCBI Gene Expression Omnibus (GEO)

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245354

This dataset is publicly available and is used exclusively for educational purposes.

---

## Repository Organization

```
rna-seq-case-study-gse245354
│
├── README.md
├── LICENSE
│
├── data
│   ├── raw
│   ├── metadata
│   └── processed
│
├── scripts
│   ├── 01_download_data.R
│   ├── 02_import_expression_matrix.R
│   ├── 03_data_cleaning.R
│   ├── 04_exploratory_analysis.R
│   ├── 05_quality_assessment.R
│   ├── 06_pca.R
│   ├── 07_heatmap.R
│   ├── 08_differential_expression.R
│   ├── 09_visualization.R
│   └── 10_export_results.R
│
├── figures
│
├── results
│
└── docs
```

---

# Course Roadmap

| Lesson | Topic |
|:-------:|-------|
| 01 | Downloading RNA-Seq data from GEO |
| 02 | Importing gene expression matrices |
| 03 | Data cleaning and organization |
| 04 | Exploratory data analysis |
| 05 | Descriptive statistics |
| 06 | Principal Component Analysis (PCA) |
| 07 | Heatmaps and hierarchical clustering |
| 08 | Differential gene expression analysis |
| 09 | Publication-quality visualizations |
| 10 | Biological interpretation of results |

---

## Learning Objectives

During this case study students will learn how to:

- Organize reproducible bioinformatics projects
- Obtain public transcriptomic datasets
- Import gene expression matrices into R
- Manipulate data using the tidyverse
- Perform exploratory data analysis
- Generate publication-quality graphics using ggplot2
- Explore sample relationships using PCA
- Perform statistical analyses
- Identify differentially expressed genes
- Interpret transcriptomic results

---

## Main R Packages

- GEOquery
- tidyverse
- ggplot2
- pheatmap
- DESeq2
- EnhancedVolcano
- patchwork

---

## Educational Philosophy

This repository was designed to bridge **R programming** and **bioinformatics** through a real-world case study.

Each script is heavily documented to explain not only *how* the code works, but also *why* each analytical step is performed, encouraging students to develop computational thinking alongside biological interpretation.

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

If you use this repository for teaching or research purposes, please cite the original GEO dataset and its associated publication.
