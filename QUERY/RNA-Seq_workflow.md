# RNA-Seq data processing and gene expression  analysis workflow

## Background
RNA seq is widely used for gene expression studies to quantify the RNA in a sample using next-generation sequencing (NGS). It is a powerful tool with many applications for gene discovery and quantification. Several tools and pipelines for RNA-seq analysis are available, and H3ABioNet provides an SOP with some recommendations for gene expression analysis in human. In the SOP, we assume differential expression is being assessed between 2 experimental conditions, i.e. a simple 1:1 comparison. The sample data is from a human genome. 

## Pipeline
The H3ABioNet have developed some [standard operating procedures (SOPs)](https://h3abionet.github.io/H3ABionet-SOPs/RNA-Seq) for the H3Africa Consortium. Although the pipeline provides a recommendation for different tools, we suggest the following for this mini-project:
- `FastQC` for quality check
- `Trimmomatic` for adaptor removal and trimming
- `HISAT2` for alignment and HTSeq’s `htseq-count`, Subread’s `featureCounts` for count generation (you are welcome to try out `salmon` and `Kallisto` for alignment-free techniques)
- `MultiQC` to collect the statistics
- use `R` for statistical analysis

The above tools are recommendations for you to start with. You are, however, encouraged to explore alternative techniques. 

## Project Task
1. Create a pipeline for RNA-seq pre-processing and gene expression analysis based on the H3ABioNet SOP
2. Convert the pipeline into Nextflow or Snakemake, and RMarkdown, making use of Conda package manager and containers
3. Document the workflow you develop on GitHub using wikis

Through the project, you need to demonstrate collaborative research skills, informative visualization, and report writing. 

