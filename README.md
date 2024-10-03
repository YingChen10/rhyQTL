# rhyQTL


[![rhyQTL](https://img.shields.io/badge/release-v1.0-brightgreen)]([https://example.com/release](https://github.com/YingChen10/rhyQTL/))
[![R](https://img.shields.io/badge/R-4.2.0-brightgreen)]([https://example.com/R](https://cran.r-project.org/))
[![bedtools](https://img.shields.io/badge/bedtools-v2.27.1-brightgreen)](https://bedtools.readthedocs.io/en/latest/)
[![samtools](https://img.shields.io/badge/samtools-v1.10-brightgreen)](https://www.htslib.org/)
[![RareAPA-Zenodo](https://img.shields.io/badge/Zenodo-brightgreen)](https://zenodo.org/records/11371829)

24-hour biological rhythms are essential to maintain physiological homeostasis. Disruption of these rhythms increases the risks of multiple diseases. The biological rhythms are known to have a genetic basis formed by core clock genes, but how individual genetic variation shapes the oscillating transcriptome and contributes to human chronophysiology and disease risk is largely unknown. Here, we mapped interactions between temporal gene expression and genotype to identify quantitative trait loci (QTLs) contributing to rhythmic gene expression. These newly identified QTLs were termed as rhythmic QTLs (rhyQTLs). Specifically, we identified rhyQTLs and rhyQTL-associated genes (rhyGenes) in 45 human tissues using data from the Genotype Tissue Expression (GTEx) Project. 

This repository contains all source code for the analyses in manuscript "[Human genetic variation determines 24-hour rhythmic gene expression and disease risk](https://www.researchsquare.com/article/rs-4790200/v1)".

## Code and data resources
This repository contains the code for data analyses described in the aforementioned rhyQTL manuscript. 
- The **Analysis** directory includes scripts for tasks such as rhyQTL mapping and enrichment analysis of rhyQTL in CREs and GWAS-tagged SNPs.Most scripts are written in R or Shell, and the required packages are specified within the code.
- Data required for the analyses can be found in the **Data** directory.

## System requirements
### Hardware and operating system requirements
A workstation or computer cluster running a POXIS system (Unix, Linux, or macOS) is required (we used the Linux distribution Ubuntu 20.04.6 LTS).
### Software and dependencies
R and dependent packages
- `R (>= 4.2.0)`
- `data.table (1.16.0)`
- `tidyverse (2.0.0)`
- `dryR (1.0.0)`
- `lmtest (0.9.40)`

External software
- `LDSC (v1.0.1)`
- `plink 1.9 beta`
- `tabix (v1.7)`
- `bedtools (>= v2.30.0)`
- `samtools (>= v1.10)`

## Authors
Ying Chen, Panpan Liu, Aniko Sabo, Dongyin Guan (Baylor College of Medicine)

## Citation
Ying Chen, Panpan Liu, Aniko Sabo, Dongyin Guan#. Human genetic variation determines 24-hour rhythmic gene expression and disease risk, 05 August 2024, Research Squar. DOI: [10.21203/rs.3.rs-4790200/v1](https://doi.org/10.21203/rs.3.rs-4790200/v1)

## Contact
If you have any comments, suggestions, questions, etc, please feel free to create a GitHub Issue.

## License
This project is covered under the MIT License.
