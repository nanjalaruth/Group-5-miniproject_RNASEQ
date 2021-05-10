[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nanjalaruth/Group-5-miniproject_RNASEQ/master)


# Group-5-miniproject_RNASEQ
As part of the EANBIT Virtual Residential Training 2020, this group will be working on a mini-project towards developing a snakemake workflow for RNA-Seq data processing and gene expression analysis. 

## Group Members
1. Ruth Nanjala (Group Lead)
1. Kakembo Fredrick Elishama
1. Eric G. Kairuki
1. Stella E. Nabirye
1. Senamile Fezile Dlamini
1. Mthande S. Mzwakhile
1. Monica Mbabazi
1. Ritah Nabunje

## A summary of the steps followed in our analysis include; 

- Pre-processing of the reads
  - Quality Check using **Fastqc** and `multiqc`
  - Trimming of poor quality bases and filtering short reads using **Trim_galore**
  - Quality check using `multiqc`
- Alignment of samples to the reference. Two approaches were used; 
  - Classical alignment using **`Hisat2`**
    - Count generation using **Subreads feature count**
      -Quality check using `multiqc`
  - Pseudo-alignment using **`Kallisto`**
- Differential Expression Analysis in R using DeSEq
- Converting the pipeline to R Markdown and Snakemake

## A summary report of the workflow is documented on the Wiki page
