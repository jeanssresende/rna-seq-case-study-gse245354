# Microarray Case Study – GSE25014

[![R](https://img.shields.io/badge/R-4.5+-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.22-019733?logo=bioconductor&logoColor=white)](https://bioconductor.org/)
[![Microarray](https://img.shields.io/badge/Microarray-Affymetrix-blueviolet)]()
[![GEO](https://img.shields.io/badge/GEO-GSE25014-success)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25014)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-In%20Development-orange)]()

---

## Overview

This repository contains a practical **microarray case study** developed to teach **R programming and transcriptomic data analysis** using the public dataset **GSE25014** from the **NCBI Gene Expression Omnibus (GEO)**.

The study investigates gene expression changes in human endothelial cells exposed to **ferric heme**, providing a real biological context for learning how to work with microarray expression data, probe-level measurements, gene annotation, exploratory analysis, visualization, and biological interpretation.

The dataset contains **24 human samples** generated using the **Affymetrix Human Genome U133 Plus 2.0 (GPL570)** platform. The experiment includes pulmonary microvascular endothelial cells (PMVECs) and pulmonary artery endothelial cells (PAECs), treated with vehicle or 5 micromolar heme. citeturn206453search0

Rather than presenting isolated R commands, this project follows a reproducible bioinformatics workflow in which programming concepts are learned through a real transcriptomic analysis.

---

## Biological Question

The original study was designed to investigate the **global gene expression response of endothelial cells to heme**, with particular interest in molecular mechanisms associated with cytoprotection and endothelial adaptation to heme-related stress. citeturn206453search0turn206453search1

For this case study, we use the dataset as a practical example to learn how to move from **microarray probes to gene-level expression** and then explore the biological differences between experimental conditions and cell types.

---

## Getting Started

Clone this repository:

```bash
git clone https://github.com/jeanssresende/microarray-case-study-gse25014.git
```

Open the project in RStudio and follow the scripts in numerical order.

The workflow will progressively cover data acquisition, import, probe annotation, probe-to-gene summarization, quality assessment, exploratory analysis, differential expression, visualization, and interpretation.

---

## Workflow

```text
GEO (GSE25014)
      |
      v
Download / Import
      |
      v
Microarray expression data
      |
      v
Probe annotation (GPL570)
      |
      v
Probe → Gene mapping
      |
      v
Probe summarization
      |
      v
Gene-level expression matrix
      |
      v
Quality assessment
      |
      v
Exploratory analysis / PCA
      |
      v
Differential expression
      |
      v
Visualization
      |
      v
Biological interpretation
```

---

## Skills Covered

- R Programming
- Data Wrangling
- Tidyverse
- Microarray Analysis
- Affymetrix Arrays
- Probe Annotation
- Probe-to-Gene Mapping
- Transcriptomics
- Exploratory Data Analysis
- Statistics
- Data Visualization
- Differential Expression Analysis
- Reproducible Research
- Bioinformatics

---

## Dataset

**Accession:** GSE25014  
**Organism:** *Homo sapiens*  
**Experiment type:** Expression profiling by array  
**Platform:** GPL570 – Affymetrix Human Genome U133 Plus 2.0 Array  
**Number of samples:** 24  
**Source:** NCBI Gene Expression Omnibus (GEO) citeturn206453search0

### Experimental groups

The dataset contains two endothelial cell models:

- Pulmonary microvascular endothelial cells (PMVECs)
- Pulmonary artery endothelial cells (PAECs)

Each cell model includes:

- Vehicle-treated samples (6 biological replicates)
- Samples treated with 5 micromolar heme (6 biological replicates)

This results in **24 samples total**. citeturn206453search0

### GEO

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25014

The data are publicly available and are used for educational and reproducibility purposes.

---

## Repository Organization

```text
microarray-case-study-gse25014
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
│   ├── 02_import_expression_data.R
│   ├── 03_probe_annotation.R
│   ├── 04_probe_to_gene.R
│   ├── 05_data_cleaning.R
│   ├── 06_quality_assessment.R
│   ├── 07_exploratory_analysis.R
│   ├── 08_pca.R
│   ├── 09_heatmap.R
│   ├── 10_differential_expression.R
│   ├── 11_visualization.R
│   └── 12_export_results.R
│
├── figures
│
├── results
│
└── docs
```

---

## Course Roadmap

| Lesson | Topic |
|:------:|-------|
| 01 | Downloading microarray data from GEO |
| 02 | Importing expression data |
| 03 | Understanding probes and microarray measurements |
| 04 | Probe annotation using GPL570 |
| 05 | Probe-to-gene mapping |
| 06 | Summarizing multiple probes per gene |
| 07 | Data cleaning and quality assessment |
| 08 | Exploratory data analysis and PCA |
| 09 | Heatmaps and sample clustering |
| 10 | Differential expression analysis |
| 11 | Publication-quality visualization |
| 12 | Biological interpretation and export |

---

## Learning Objectives

During this case study students will learn how to:

- Organize a reproducible bioinformatics project
- Obtain public transcriptomic datasets from GEO
- Understand the structure of Affymetrix microarray data
- Work with probe-level expression measurements
- Annotate probes using platform information
- Convert probe identifiers into gene identifiers
- Resolve situations where multiple probes map to the same gene
- Build a gene-level expression matrix
- Perform quality assessment and exploratory analysis
- Explore sample relationships using PCA and clustering
- Perform differential expression analysis
- Generate publication-quality graphics
- Interpret transcriptomic results in a biological context

---

## Main R Packages

- GEOquery
- Bioconductor
- tidyverse
- limma
- annotate
- hgu133plus2.db
- ggplot2
- pheatmap
- EnhancedVolcano
- patchwork

---

## Educational Philosophy

This repository was designed to bridge **R programming**, **microarray analysis**, and **bioinformatics** through a real-world case study.

Each script is documented to explain not only *how* the code works, but also *why* each analytical step is performed. Special attention is given to the transition from **probe-level measurements to gene-level expression**, an essential concept for students working with legacy microarray datasets.

The goal is to develop computational thinking together with biological interpretation and reproducible research practices.

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

Please cite the original **GSE25014** GEO record and its associated publication when using this repository for teaching or research purposes. The GEO record reports the study *Gene expression data of endothelium exposed to heme*, and links the dataset to the publication by Ghosh et al. (2011). citeturn206453search0turn206453search1
