# Worm Perturb-Seq (WPS) FASTQ Pipeline Evaluation

## Overview
This task is based on the dataset from the Worm Perturb-Seq project (2024), focusing on high-throughput single-cell transcriptomics in *C. elegans* using barcoded CRISPR screening.

## Objective
The goal is to evaluate the ability of a system to reproduce a FASTQ-based pipeline using a small sample of data (SRR25848640). The questions are designed around different steps of the pipeline from quality control to gene quantification.

## Dataset
- **BioProject**: PRJNA1066733
- **GEO Accession**: GSE253847
- **Sample**: SRR25848640 (subsampled)
- **Data Size**: <100MB (subsampled FASTQ)

## Pipeline Steps Covered
- Quality Control
- Alignment (to *C. elegans* genome)
- UMI Deduplication
- Cell Filtering
- Gene Quantification

## How to Use
1. Run the provided subsampled FASTQ files through the Worm Perturb-Seq pipeline available on GitHub.
2. Answer the 5 objective questions provided in `questions.yaml`.
3. Validate your results against the provided ground truth in `answers.yaml`.

## References
- **Paper**: "Combinatorial genetics reveals interactions between splicing regulators in C. elegans" (2024)
- **DOI**: [10.1101/2024.03.08.583410](https://doi.org/10.1101/2024.03.08.583410)
- **GitHub**: [https://github.com/AllenInstitute/worm-perturb-seq](https://github.com/AllenInstitute/worm-perturb-seq)
